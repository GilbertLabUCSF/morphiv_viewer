"""EBs-first gene profile for the MORPHIC TF perturbation screen."""

from pathlib import Path
import sys

import pandas as pd
import streamlit as st


st.set_page_config(
    page_title="Gene Profile · MORPHIC EBs",
    page_icon="logo_transparent.ico",
    layout="wide",
)

sys.path.insert(0, str(Path(__file__).parent.parent))

from src.components import render_sidebar
from src.config import (
    DEG_TABLES_EBS,
    get_antibody_validation_path,
    get_deg_volcano_path,
    get_dose_response_panel_path,
    get_edist_figure_path,
    get_edist_histogram_path,
    get_marker_dotplot_path,
    get_marker_tpm_path,
    get_pathway_heatmap_path,
    get_spider_plot_path,
    get_tf_similarity_path,
    get_umap_highlight_path,
    get_viability_figure_path,
)
from src.data_loader import (
    get_available_lineage_de,
    get_gene_list,
    get_perturbations_for_gene,
    load_compositional,
    load_deg_table,
    load_dose_response,
    load_knockdown_efficiency,
    load_lineage_analysis,
    load_lineage_de,
    load_pathway_enrichment_for_perturbation,
    load_probability_shifts,
    load_tf_clusters,
    load_timecourse_expression,
    load_transcriptome_edist,
    load_viability,
)
from src.data_model import lineage_label, perturbation_profile, timepoint_label
from src.plots import (
    plot_compositional_bars,
    plot_deg_volcano,
    plot_dose_response_scatter,
    plot_edist_rank,
    plot_lineage_de_volcano,
    plot_lineage_effect_bars,
    plot_pathway_enrichment_bars,
    plot_probability_shift_radar,
    plot_timecourse_expression,
    plot_viability_rank,
)
from src.styles import inject_css


inject_css()

try:
    genes = get_gene_list()
    lineage_df = load_lineage_analysis()
except FileNotFoundError:
    st.error("The EBs screen export is unavailable.")
    render_sidebar()
    st.stop()

query_gene = st.query_params.get("gene", "")
if isinstance(query_gene, list):
    query_gene = query_gene[0] if query_gene else ""

if query_gene in genes and st.session_state.get("_route_gene") != query_gene:
    st.session_state["gene_detail_selector"] = query_gene
    st.session_state["_route_gene"] = query_gene

st.markdown('<p class="eyebrow">EBs GENE PROFILE</p>', unsafe_allow_html=True)
selector_columns = st.columns([1.5, 1])
with selector_columns[0]:
    gene = st.selectbox(
        "Gene",
        options=[""] + genes,
        format_func=lambda value: "Choose a gene…" if not value else value,
        key="gene_detail_selector",
    )

if not gene:
    st.title("Gene profile")
    st.info("Choose a gene to inspect its EBs phenotype and timecourse context.")
    render_sidebar()
    st.stop()

st.session_state["selected_gene"] = gene
if st.query_params.get("gene") != gene:
    st.query_params["gene"] = gene

gene_rows = lineage_df[lineage_df["gene"] == gene].copy()
perturbations = get_perturbations_for_gene(gene)["perturbations"]
if not perturbations:
    st.error("No EBs perturbation profile is available for this gene.")
    render_sidebar(current_gene=gene)
    st.stop()

query_perturbation = st.query_params.get("pert", "")
if isinstance(query_perturbation, list):
    query_perturbation = query_perturbation[0] if query_perturbation else ""
default_perturbation = (
    query_perturbation if query_perturbation in perturbations else perturbations[0]
)

with selector_columns[1]:
    selected_perturbation = st.selectbox(
        "Perturbation",
        options=perturbations,
        index=perturbations.index(default_perturbation),
    )

st.markdown(
    f'<h1 class="profile-heading"><span>{gene}</span>'
    '<small>Embryoid-body CRISPRi phenotype · '
    "developmental expression context</small></h1>",
    unsafe_allow_html=True,
)

if st.query_params.get("pert") != selected_perturbation:
    st.query_params["pert"] = selected_perturbation

selected_rows = gene_rows[gene_rows["perturbation"] == selected_perturbation]
profile = perturbation_profile(gene_rows, selected_perturbation)

try:
    probability_shifts = load_probability_shifts()
    selected_shifts = probability_shifts[
        probability_shifts["perturbation"] == selected_perturbation
    ].copy()
except FileNotFoundError:
    probability_shifts = pd.DataFrame()
    selected_shifts = pd.DataFrame()

paper_shifts = selected_shifts[
    ~selected_shifts.get("lineage", pd.Series(dtype=str))
    .astype(str)
    .str.lower()
    .isin({"unknown", "unassigned"})
].copy()
if not paper_shifts.empty:
    paper_shifts["log2fc"] = pd.to_numeric(
        paper_shifts["log2fc"], errors="coerce"
    )
    paper_shifts = paper_shifts.dropna(subset=["log2fc"])

if not paper_shifts.empty:
    strongest_knn_row = paper_shifts.loc[paper_shifts["log2fc"].abs().idxmax()]
    strongest_knn_lfc = float(strongest_knn_row["log2fc"])
    strongest_knn_lineage = lineage_label(str(strongest_knn_row["lineage"]))
    if "significant" in paper_shifts:
        significant_mask = paper_shifts["significant"].fillna(False)
        if significant_mask.dtype != bool:
            significant_mask = (
                significant_mask.astype(str).str.lower().isin({"true", "1", "yes"})
            )
    elif "q_value_boot" in paper_shifts:
        significant_mask = paper_shifts["q_value_boot"].le(0.05).fillna(False)
    else:
        significant_mask = pd.Series(False, index=paper_shifts.index)
    significant_knn_shifts = paper_shifts[significant_mask].sort_values(
        "log2fc", key=lambda values: values.abs(), ascending=False
    )
else:
    strongest_knn_lfc = None
    strongest_knn_lineage = ""
    significant_knn_shifts = pd.DataFrame()

try:
    knockdown_df = load_knockdown_efficiency()
    knockdown_row = knockdown_df[
        knockdown_df["perturbation"] == selected_perturbation
    ]
except FileNotFoundError:
    knockdown_row = pd.DataFrame()

try:
    viability_df = load_viability()
    viability_row = viability_df[viability_df["gene"] == gene]
except FileNotFoundError:
    viability_df = pd.DataFrame()
    viability_row = pd.DataFrame()

knockdown_value = (
    float(knockdown_row.iloc[0]["knockdown_pct"]) if not knockdown_row.empty else None
)
viability_value = (
    float(viability_row.iloc[0]["median_lfc"]) if not viability_row.empty else None
)

profile_metrics = [
    ("Cells", f"{profile.get('n_cells', 0):,}", ""),
    (
        "Strongest KNN shift",
        f"{strongest_knn_lfc:+.2f}" if strongest_knn_lfc is not None else "N/A",
        f"{strongest_knn_lineage} · log2FC" if strongest_knn_lineage else "",
    ),
    ("Significant KNN shifts", str(len(significant_knn_shifts)), "q ≤ 0.05"),
    (
        "Target knockdown",
        f"{knockdown_value:.1f}%" if knockdown_value is not None else "N/A",
        "",
    ),
    (
        "Viability LFC",
        f"{viability_value:.2f}" if viability_value is not None else "N/A",
        "",
    ),
]
metric_markup = "".join(
    f'<div class="metric-card"><div class="profile-metric-label">{label}</div>'
    f'<div class="profile-metric-value">{value}</div>'
    f'<div class="profile-metric-note">{note}</div></div>'
    for label, value, note in profile_metrics
)
st.markdown(
    f'<div class="profile-metric-grid">{metric_markup}</div>',
    unsafe_allow_html=True,
)

if not significant_knn_shifts.empty:
    strongest_effect = significant_knn_shifts.iloc[0]
    effect = float(strongest_effect["log2fc"])
    strongest_q = pd.to_numeric(
        pd.Series([strongest_effect.get("q_value_boot")]), errors="coerce"
    ).iloc[0]
    strongest_q_label = f"{strongest_q:.2g}" if pd.notna(strongest_q) else "N/A"
    direction = "increased" if effect > 0 else "decreased"
    additional = len(significant_knn_shifts) - 1
    additional_text = (
        f"; {additional} additional lineage{'s' if additional != 1 else ''} also pass q ≤ 0.05"
        if additional
        else ""
    )
    st.markdown(
        f'<div class="evidence-callout"><strong>Primary readout:</strong> '
        f"{lineage_label(str(strongest_effect['lineage']))} mean KNN probability "
        f"is {direction} (NTC-centered log2FC {effect:+.2f}, "
        f"q = {strongest_q_label}){additional_text}.</div>",
        unsafe_allow_html=True,
    )
elif not paper_shifts.empty:
    st.markdown(
        '<div class="evidence-callout"><strong>Primary readout:</strong> '
        "no KNN lineage-probability shift passes q ≤ 0.05 for this "
        "perturbation.</div>",
        unsafe_allow_html=True,
    )
else:
    st.markdown(
        '<div class="evidence-callout"><strong>Primary readout:</strong> '
        "the KNN probability-shift table is unavailable for this perturbation.</div>",
        unsafe_allow_html=True,
    )
st.divider()

view = st.segmented_control(
    "Profile view",
    options=[
        "Phenotype",
        "Expression over time",
        "Molecular readouts",
        "Data & methods",
    ],
    default="Phenotype",
    key="profile_view",
    label_visibility="collapsed",
)


if view == "Phenotype":
    st.subheader("Lineage phenotype")
    spider_path = get_spider_plot_path(selected_perturbation, "EBs")
    umap_path = get_umap_highlight_path(selected_perturbation, "EBs")

    if not selected_shifts.empty:
        if umap_path.exists():
            spider_column, highlight_column = st.columns(2)
            with spider_column:
                st.plotly_chart(
                    plot_probability_shift_radar(
                        probability_shifts, selected_perturbation
                    ),
                    width="stretch",
                    key=f"interactive_knn_lfc_{selected_perturbation}",
                    config={
                        "displayModeBar": True,
                        "displaylogo": False,
                        "responsive": True,
                    },
                )
            with highlight_column:
                st.image(
                    str(umap_path),
                    caption=(
                        f"Pipeline-generated EBs embedding highlight for "
                        f"{selected_perturbation}"
                    ),
                    width="stretch",
                )
            st.caption(
                "Left: interactive rendering of the paper's NTC-centered log2 fold "
                "change in mean KNN lineage probability. Right: the matching "
                "pipeline-generated perturbation-cell highlight."
            )
        else:
            _, spider_column, _ = st.columns([0.7, 2, 0.7])
            with spider_column:
                st.plotly_chart(
                    plot_probability_shift_radar(
                        probability_shifts, selected_perturbation
                    ),
                    width="stretch",
                    key=f"interactive_knn_lfc_{selected_perturbation}",
                    config={
                        "displayModeBar": True,
                        "displaylogo": False,
                        "responsive": True,
                    },
                )
            st.caption(
                "Interactive rendering of the paper's NTC-centered log2 fold change "
                "in mean KNN lineage probability. Hover over a lineage marker for "
                "its confidence interval, q-value, probabilities, and cell count."
            )
    elif spider_path.exists():
        st.warning(
            "The KNN estimates are unavailable; showing the pipeline-rendered spider "
            "plot as a fallback."
        )
        _, spider_column, _ = st.columns([0.7, 2, 0.7])
        with spider_column:
            st.image(
                str(spider_path),
                caption=f"{selected_perturbation}: pipeline-generated spider plot",
                width="stretch",
            )
        st.caption(
            "Static fallback using the manuscript pipeline's NTC-centered KNN log2FC."
        )
    else:
        st.info("A spider plot is not available for this perturbation.")

    if not selected_shifts.empty:
        shift_columns = [
            "lineage",
            "log2fc",
            "ci_low",
            "ci_high",
            "q_value_boot",
            "significant",
            "mean_prob_pert",
            "mean_prob_ntc",
            "n_cells",
        ]
        shift_table = selected_shifts[
            [column for column in shift_columns if column in selected_shifts]
        ].copy()
        shift_table["lineage"] = shift_table["lineage"].map(lineage_label)
        shift_display = shift_table.copy()
        for probability_column in ["mean_prob_pert", "mean_prob_ntc"]:
            if probability_column in shift_display:
                shift_display[probability_column] = (
                    pd.to_numeric(shift_display[probability_column], errors="coerce")
                    * 100
                )
        shift_display = shift_display.rename(
            columns={
                "lineage": "Lineage",
                "log2fc": "NTC-centered log2FC",
                "ci_low": "CI low",
                "ci_high": "CI high",
                "q_value_boot": "q-value",
                "significant": "Significant",
                "mean_prob_pert": "Perturbed mean probability",
                "mean_prob_ntc": "NTC mean probability",
                "n_cells": "Cells",
            }
        )
        with st.expander("View and download the underlying KNN estimates"):
            st.dataframe(
                shift_display,
                column_config={
                    "NTC-centered log2FC": st.column_config.NumberColumn(format="%.2f"),
                    "CI low": st.column_config.NumberColumn(format="%.2f"),
                    "CI high": st.column_config.NumberColumn(format="%.2f"),
                    "q-value": st.column_config.NumberColumn(format="%.3g"),
                    "Perturbed mean probability": st.column_config.NumberColumn(
                        format="%.2f%%"
                    ),
                    "NTC mean probability": st.column_config.NumberColumn(
                        format="%.2f%%"
                    ),
                    "Cells": st.column_config.NumberColumn(format="%d"),
                },
                hide_index=True,
                width="stretch",
            )
            st.download_button(
                "Download KNN log2FC estimates",
                data=shift_table.to_csv(index=False),
                file_name=f"{selected_perturbation}_knn_log2fc.csv",
                mime="text/csv",
            )

    with st.expander("Supplementary lineage-score view (Glass's Δ)"):
        effect_figure = plot_lineage_effect_bars(gene_rows, selected_perturbation)
        st.plotly_chart(effect_figure, width="stretch")
        st.caption(
            "This is a separate lineage-score effect size, not the KNN probability "
            "log2FC used by the manuscript spider plot. Faded bars do not pass q ≤ 0.05."
        )

        evidence = selected_rows[
            [
                "lineage",
                "observed_glass_delta",
                "q_value",
                "effect_size_interpretation",
                "n_cells",
            ]
        ].copy()
        evidence["lineage"] = evidence["lineage"].map(lineage_label)
        evidence = evidence.rename(
            columns={
                "lineage": "Lineage",
                "observed_glass_delta": "Glass's Δ",
                "q_value": "q-value",
                "effect_size_interpretation": "Effect class",
                "n_cells": "Cells",
            }
        ).sort_values("Glass's Δ", key=lambda values: values.abs(), ascending=False)
        st.dataframe(
            evidence,
            column_config={
                "Glass's Δ": st.column_config.NumberColumn(format="%.3f"),
                "q-value": st.column_config.NumberColumn(format="%.3g"),
                "Cells": st.column_config.NumberColumn(format="%d"),
            },
            hide_index=True,
            width="stretch",
        )
        st.download_button(
            "Download lineage estimates",
            data=evidence.to_csv(index=False),
            file_name=f"{selected_perturbation}_ebs_lineage_effects.csv",
            mime="text/csv",
        )

    st.subheader("Supporting phenotype context")
    viability_source = get_viability_figure_path("EBs")
    if viability_source.exists():
        with st.expander("Viability / fitness source figure"):
            st.image(
                str(viability_source),
                caption="Pipeline-generated EBs viability volcano",
                width="stretch",
            )
            if not viability_df.empty:
                st.plotly_chart(
                    plot_viability_rank(viability_df, gene, "EBs"),
                    width="stretch",
                )
                st.caption(
                    f"Interactive locator for {gene}; the volcano above is the "
                    "authoritative pipeline rendering."
                )
    elif not viability_df.empty:
        st.plotly_chart(
            plot_viability_rank(viability_df, gene, "EBs"), width="stretch"
        )
        st.caption(
            "The pipeline viability figure is unavailable; this interactive rank "
            "is derived from the exported viability table."
        )

    try:
        compositional_df = load_compositional()
        gene_compositional = compositional_df[compositional_df["gene"] == gene]
    except FileNotFoundError:
        gene_compositional = pd.DataFrame()
    if not gene_compositional.empty:
        with st.expander("Significant compositional probability shifts"):
            st.plotly_chart(
                plot_compositional_bars(gene_compositional, gene),
                width="stretch",
            )


elif view == "Expression over time":
    st.subheader("Baseline expression through differentiation")
    st.markdown(
        "This timecourse is an independent, unperturbed reference. It shows when the target "
        "is normally expressed; it does not measure the consequence of its knockdown."
    )
    try:
        timecourse_df = load_timecourse_expression()
        gene_timecourse = timecourse_df[timecourse_df["gene"] == gene].copy()
    except FileNotFoundError:
        gene_timecourse = pd.DataFrame()

    if gene_timecourse.empty:
        st.info(f"{gene} is not present in the current timecourse expression export.")
    else:
        peak = gene_timecourse.loc[gene_timecourse["mean_expression"].idxmax()]
        timecourse_metrics = st.columns(3)
        timecourse_metrics[0].metric("Peak timepoint", timepoint_label(peak["day"]))
        timecourse_metrics[1].metric("Peak mean expression", f"{peak['mean_expression']:.2f}")
        detected_column = "pct_detected" if "pct_detected" in peak else "pct_expressing"
        timecourse_metrics[2].metric("Cells detected at peak", f"{peak[detected_column]:.1f}%")

        st.plotly_chart(
            plot_timecourse_expression(timecourse_df, gene),
            width="stretch",
        )
        st.caption(
            "Mean expression is log1p-normalized counts per 10,000. Detection is the "
            "percentage of cells with at least one observed count."
        )
        with st.expander("View and download timecourse values"):
            st.dataframe(gene_timecourse, hide_index=True, width="stretch")
            st.download_button(
                "Download timecourse values",
                data=gene_timecourse.to_csv(index=False),
                file_name=f"{gene}_timecourse_expression.csv",
                mime="text/csv",
            )

    try:
        timecourse_url = st.secrets.get("cellxgene", {}).get("timecourse_url", "")
    except Exception:
        timecourse_url = ""
    if timecourse_url:
        st.link_button(
            "Explore the full timecourse in CellxGene",
            timecourse_url,
            icon=":material/open_in_new:",
        )


elif view == "Molecular readouts":
    st.subheader("Supporting molecular readouts")
    analysis = st.selectbox(
        "Analysis",
        options=[
            "Differential expression",
            "Lineage-specific differential expression",
            "Pathway enrichment",
            "Transcriptome shift and dose response",
            "Marker validation",
            "TF similarity",
        ],
        key="molecular_analysis",
    )

    if analysis == "Differential expression":
        deg_df = load_deg_table(selected_perturbation)
        volcano_path = get_deg_volcano_path(selected_perturbation, "EBs")
        if deg_df is not None:
            st.plotly_chart(
                plot_deg_volcano(deg_df, f"{selected_perturbation} · EBs"),
                width="stretch",
                key=f"interactive_deg_volcano_{selected_perturbation}",
                config={
                    "displayModeBar": True,
                    "displaylogo": False,
                    "responsive": True,
                },
            )
            st.caption(
                "Interactive rendering from the same exported DEG table and significance "
                "thresholds as the pipeline figure. Hover over a point for rounded gene-level "
                "statistics."
            )
        elif volcano_path.exists():
            st.warning(
                "The DEG table is unavailable; showing the pipeline volcano as a fallback."
            )
            st.image(
                str(volcano_path),
                caption=(
                    f"{selected_perturbation}: pipeline-generated EBs DEG volcano"
                ),
                width="stretch",
            )
        else:
            st.info("Differential-expression results are not available for this perturbation.")

        if deg_df is not None:
            significant_degs = (
                deg_df[deg_df["padj"] <= 0.05].sort_values("padj")
                if "padj" in deg_df else deg_df
            )
            deg_number_formats = {
                "baseMean": "%.1f",
                "log2FoldChange": "%.2f",
                "lfcSE": "%.2f",
                "stat": "%.2f",
                "pvalue": "%.3g",
                "padj": "%.3g",
            }
            deg_column_config = {
                column: st.column_config.NumberColumn(format=number_format)
                for column, number_format in deg_number_formats.items()
                if column in significant_degs
            }
            with st.expander("View the top significant genes"):
                st.dataframe(
                    significant_degs.head(50),
                    column_config=deg_column_config,
                    hide_index=True,
                    width="stretch",
                )
            st.download_button(
                "Download DEG table",
                data=deg_df.to_csv(index=False),
                file_name=f"{selected_perturbation}_ebs_degs.csv",
                mime="text/csv",
            )

    elif analysis == "Lineage-specific differential expression":
        available_lineages = get_available_lineage_de(selected_perturbation)
        if not available_lineages:
            st.info("No lineage-specific DE result is available for this perturbation.")
        else:
            selected_lineage = st.selectbox(
                "Lineage",
                options=available_lineages,
                format_func=lineage_label,
                key="lineage_de_selector",
            )
            lineage_de = load_lineage_de(selected_perturbation, selected_lineage)
            if lineage_de is not None:
                st.plotly_chart(
                    plot_lineage_de_volcano(
                        lineage_de, selected_perturbation, lineage_label(selected_lineage)
                    ),
                    width="stretch",
                )
                st.download_button(
                    "Download lineage DE table",
                    data=lineage_de.to_csv(index=False),
                    file_name=f"{selected_perturbation}_{selected_lineage}_ebs_de.csv",
                    mime="text/csv",
                )

    elif analysis == "Pathway enrichment":
        pathway_source = get_pathway_heatmap_path("EBs")
        if pathway_source.exists():
            st.image(
                str(pathway_source),
                caption="Pipeline-generated EBs pathway-enrichment overview",
                width="stretch",
            )
            st.caption("Authoritative overview from the analysis pipeline.")

        with st.spinner("Loading enrichment results for this perturbation…"):
            enrichment = load_pathway_enrichment_for_perturbation(selected_perturbation)
        if enrichment is None or enrichment.empty:
            st.info("No significant pathway enrichment is available for this perturbation.")
        else:
            libraries = sorted(enrichment["gene_set_library"].dropna().unique())
            selected_library = st.selectbox("Gene-set library", ["All"] + libraries)
            displayed_enrichment = (
                enrichment if selected_library == "All"
                else enrichment[enrichment["gene_set_library"] == selected_library]
            )
            with st.expander(
                f"Per-perturbation browser view for {selected_perturbation}",
                expanded=not pathway_source.exists(),
            ):
                st.plotly_chart(
                    plot_pathway_enrichment_bars(
                        displayed_enrichment, selected_perturbation
                    ),
                    width="stretch",
                )
                st.caption(
                    "Interactive subset derived from the exported pathway table; it is "
                    "supplementary to the pipeline overview above."
                )
                st.dataframe(displayed_enrichment, hide_index=True, width="stretch")
            st.download_button(
                "Download enrichment results",
                data=displayed_enrichment.to_csv(index=False),
                file_name=f"{selected_perturbation}_ebs_pathways.csv",
                mime="text/csv",
            )

    elif analysis == "Transcriptome shift and dose response":
        edist_source = get_edist_figure_path("EBs")
        edist_histogram = get_edist_histogram_path("EBs")
        dose_source = get_dose_response_panel_path("EBs")
        source_figures = [
            (edist_source, "Pipeline E-distance top-20 overview"),
            (edist_histogram, "Pipeline E-distance distribution"),
            (dose_source, "Pipeline dose-response panel"),
        ]
        available_sources = [item for item in source_figures if item[0].exists()]
        if available_sources:
            source_columns = st.columns(len(available_sources))
            for column, (source_path, caption) in zip(
                source_columns, available_sources
            ):
                with column:
                    st.image(str(source_path), caption=caption, width="stretch")

        edist = load_transcriptome_edist()
        dose_response = load_dose_response()
        if edist is not None or dose_response is not None:
            with st.expander(
                f"Interactive locators for {gene}", expanded=not available_sources
            ):
                result_columns = st.columns(2)
                with result_columns[0]:
                    if edist is not None:
                        st.plotly_chart(
                            plot_edist_rank(edist, gene, "EBs"), width="stretch"
                        )
                    else:
                        st.info("Transcriptome E-distance is unavailable.")
                with result_columns[1]:
                    if dose_response is not None:
                        st.plotly_chart(
                            plot_dose_response_scatter(
                                dose_response, gene, "EBs"
                            ),
                            width="stretch",
                        )
                    else:
                        st.info("Dose-response results are unavailable.")
                st.caption(
                    "These interactive locators are derived browser views; the panels "
                    "above are the pipeline-generated source figures."
                )

    elif analysis == "Marker validation":
        marker_dotplot = get_marker_dotplot_path(selected_perturbation, "EBs")
        marker_cpm = get_marker_tpm_path(selected_perturbation, "EBs")
        antibody = get_antibody_validation_path(selected_perturbation, "EBs")

        marker_images = [
            (marker_dotplot, "Marker-gene dot plot"),
            (marker_cpm, "Perturbed versus NTC marker expression · CPM"),
        ]
        shown = False
        for image_path, caption in marker_images:
            if image_path.exists():
                st.image(str(image_path), caption=caption, width="stretch")
                shown = True

        if antibody.exists():
            st.markdown("#### Antibody-marker panel")
            st.image(
                str(antibody),
                caption="Perturbed versus NTC expression for antibody-detectable markers · CPM",
                width="stretch",
            )
            st.caption(
                "This panel restricts the expression comparison to markers with antibodies "
                "selected for experimental validation."
            )
            shown = True

        if not shown:
            st.info("No marker-validation figure is available for this perturbation.")

    elif analysis == "TF similarity":
        clusters = load_tf_clusters()
        gene_cluster = (
            clusters[clusters["perturbation"] == selected_perturbation]
            if clusters is not None else pd.DataFrame()
        )
        if gene_cluster.empty:
            st.info("No TF-similarity cluster is available for this perturbation.")
        else:
            st.metric("Transcriptomic signature cluster", int(gene_cluster.iloc[0]["cluster"]))
            similarity_figure = get_tf_similarity_path("EBs")
            if similarity_figure.exists():
                st.image(
                    str(similarity_figure),
                    caption="EBs TF-similarity overview",
                    width="stretch",
                )


elif view == "Data & methods":
    st.subheader("What is available for this perturbation?")
    try:
        probability_perturbations = set(load_probability_shifts()["perturbation"])
    except FileNotFoundError:
        probability_perturbations = set()
    availability = {
        "Pipeline spider figure": get_spider_plot_path(
            selected_perturbation, "EBs"
        ).exists(),
        "KNN probability-shift table": selected_perturbation
        in probability_perturbations,
        "UMAP highlight": get_umap_highlight_path(selected_perturbation, "EBs").exists(),
        "Pipeline DEG volcano": get_deg_volcano_path(
            selected_perturbation, "EBs"
        ).exists(),
        "Differential-expression table": (
            DEG_TABLES_EBS / f"{selected_perturbation}_deg.csv"
        ).exists(),
        "Lineage-specific DE": bool(
            get_available_lineage_de(selected_perturbation)
        ),
        "Marker expression (CPM)": get_marker_dotplot_path(
            selected_perturbation, "EBs"
        ).exists() or get_marker_tpm_path(selected_perturbation, "EBs").exists(),
        "Antibody-marker panel (CPM)": get_antibody_validation_path(
            selected_perturbation, "EBs"
        ).exists(),
    }
    availability_table = pd.DataFrame(
        {
            "Analysis": availability.keys(),
            "Status": ["Available" if value else "Not available" for value in availability.values()],
        }
    )
    st.dataframe(availability_table, hide_index=True, width="stretch")

    st.markdown(
        """
        ### Reading the profile

        - **Glass's Δ** is the standardized shift in a lineage score relative to non-targeting controls. Sign indicates direction; magnitude indicates effect size.
        - **q-value** is the multiple-testing-adjusted significance value. This portal marks q ≤ 0.05 as significant.
        - **Spider plots** use the pipeline's NTC-centered log2 fold change in mean KNN lineage probability. The primary interactive view uses the same exported values, significance calls, thresholds, and scale as the manuscript pipeline.
        - **Figure policy:** the spider and DEG volcano are interactive primary views because they reproduce the pipeline data and rules directly. Their static pipeline assets are fallbacks; other manuscript figures remain direct source images unless an equivalent interactive rendering exists.
        - **Viability LFC** describes representation of a perturbed gene relative to expectation. Strong depletion can confound downstream phenotype interpretation.
        - **Timecourse expression** comes from a separate, unperturbed differentiation experiment and supplies developmental context only.
        """
    )
    st.caption(
        "The portal intentionally focuses on the EBs screen in this release. Analyses are "
        "shown only when their current export exists."
    )


render_sidebar(current_gene=gene)
