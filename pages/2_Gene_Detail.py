"""
Gene Detail Page - Comprehensive view of a single TF perturbation
"""

import streamlit as st
import pandas as pd

st.set_page_config(
    page_title="Gene Detail - MORPHIC Portal",
    page_icon="logo_transparent.ico",
    layout="wide",
)

# Add src to path
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent))

from src.styles import inject_css
from src.components import render_sidebar
inject_css()

from src.data_loader import (
    get_gene_list,
    get_perturbations_for_gene,
    load_lineage_analysis,
    load_knockdown_efficiency,
    load_viability,
    load_compositional,
    load_deg_table,
    load_timecourse_expression,
    load_dose_response,
    load_pathway_enrichment,
    load_transcriptome_edist,
    load_tf_clusters,
    load_lineage_de,
    get_available_lineage_de,
)
from src.plots import (
    plot_lineage_heatmap,
    plot_knockdown_comparison,
    plot_viability_rank,
    plot_deg_volcano,
    plot_timecourse_expression,
    plot_compositional_bars,
    plot_edist_rank,
    plot_dose_response_scatter,
    plot_pathway_enrichment_bars,
    plot_lineage_de_volcano,
)
from src.config import (
    get_umap_highlight_path,
    get_marker_dotplot_path,
    get_spider_plot_path,
    get_marker_tpm_path,
    get_antibody_validation_path,
    get_tf_similarity_path,
    CONDITIONS,
)
from src.gene_summary import render_gene_summary_section

st.title("Gene Detail")

# Gene selector
try:
    genes = get_gene_list()
except FileNotFoundError as e:
    st.error(f"Cannot load gene list: {e}")
    st.stop()

# Accept gene from query params (?gene=SOX2) or session
query_params = st.query_params if hasattr(st, "query_params") else st.experimental_get_query_params()
qp_gene = query_params.get("gene", "")
if isinstance(qp_gene, list):
    qp_gene = qp_gene[0] if qp_gene else ""

default_gene = qp_gene or st.session_state.get("selected_gene", "")

if default_gene and default_gene in genes:
    default_idx = genes.index(default_gene) + 1  # +1 because of empty option
else:
    default_idx = 0

gene = st.selectbox(
    "Select a gene:",
    options=[""] + genes,
    index=default_idx,
    format_func=lambda x: "Choose a gene..." if x == "" else x,
    key="gene_detail_selector",
)

if not gene:
    st.info("Select a gene from the dropdown above to view details.")
    st.stop()

# Store selection
st.session_state["selected_gene"] = gene

# Header with gene name
st.header(gene)

# Helper function for knockdown category
def get_knockdown_category(pct):
    if pct >= 90:
        return "Excellent"
    elif pct >= 70:
        return "Good"
    elif pct >= 30:
        return "Moderate"
    else:
        return "Poor"

# Load all metrics data upfront
glass_delta_max = None
via_values = {"EBs": None, "iPSC": None}
kd_values = {"iPSC": None, "EBs": None}

try:
    _lineage_df = load_lineage_analysis()
    _gene_lineage = _lineage_df[_lineage_df["gene"] == gene]
    if not _gene_lineage.empty and "observed_glass_delta" in _gene_lineage.columns:
        glass_delta_max = _gene_lineage["observed_glass_delta"].abs().max()
except FileNotFoundError:
    pass

try:
    _via_df = load_viability()
    _via_gene = _via_df[_via_df["gene"] == gene]
    for cond in ["EBs", "iPSC"]:
        row = _via_gene[_via_gene["condition"] == cond]
        if not row.empty:
            via_values[cond] = row["median_lfc"].iloc[0]
except FileNotFoundError:
    pass

try:
    kd_df = load_knockdown_efficiency()
    gene_kd = kd_df[kd_df["target_gene"] == gene]
    for cond in ["iPSC", "EBs"]:
        cond_data = gene_kd[gene_kd["condition"] == cond]
        if len(cond_data) > 0:
            kd_values[cond] = cond_data.iloc[0]["knockdown_pct"]
except Exception:
    pass

# Display all metrics in aligned 5-column layout
col1, col2, col3, col4, col5 = st.columns(5)

with col1:
    st.metric("Glass's \u0394 (max)", f"{glass_delta_max:.2f}" if glass_delta_max is not None else "N/A")

with col2:
    val = via_values["iPSC"]
    st.metric("Viability (iPSC)", f"{val:.2f}" if val is not None else "N/A")

with col3:
    val = via_values["EBs"]
    st.metric("Viability (EBs)", f"{val:.2f}" if val is not None else "N/A")

with col4:
    val = kd_values["iPSC"]
    if val is not None:
        st.metric("Knockdown (iPSC)", f"{val:.1f}%", delta=get_knockdown_category(val), delta_color="off")
    else:
        st.metric("Knockdown (iPSC)", "N/A")

with col5:
    val = kd_values["EBs"]
    if val is not None:
        st.metric("Knockdown (EBs)", f"{val:.1f}%", delta=get_knockdown_category(val), delta_color="off")
    else:
        st.metric("Knockdown (EBs)", "N/A")

st.divider()

# Tabs for different views
tab_overview, tab_perturbation, tab_lineage_de, tab_pathway, tab_markers, tab_timecourse, tab_cellxgene = st.tabs([
    "Overview",
    "Perturbation Effects",
    "Lineage DE",
    "Pathway Enrichment",
    "Marker Expression",
    "Timecourse",
    "CellxGene",
])

# ============ OVERVIEW TAB ============
with tab_overview:
    render_gene_summary_section(gene)


# ============ PERTURBATION TAB ============
with tab_perturbation:
    # Get perturbation IDs for this gene (used in multiple sections)
    try:
        pert_info = get_perturbations_for_gene(gene)
        perturbations = pert_info["perturbations"]
    except Exception:
        perturbations = []

    # UMAP Highlight section
    st.subheader("Cell Distribution (UMAP)")

    if perturbations:
        umap_pert = st.selectbox(
            "Select perturbation:",
            options=perturbations,
            key="umap_pert_select",
        )

        umap_col1, umap_col2 = st.columns(2)

        with umap_col1:
            st.markdown("**EBs**")
            umap_path_ebs = get_umap_highlight_path(umap_pert, "EBs")
            if umap_path_ebs.exists():
                st.image(str(umap_path_ebs), use_container_width=True)
            else:
                st.info(f"No UMAP highlight for {umap_pert} in EBs")

        with umap_col2:
            st.markdown("**iPSC**")
            umap_path_ipsc = get_umap_highlight_path(umap_pert, "iPSC")
            if umap_path_ipsc.exists():
                st.image(str(umap_path_ipsc), use_container_width=True)
            else:
                st.info(f"No UMAP highlight for {umap_pert} in iPSC")
    else:
        st.warning(f"No perturbations found for {gene}")

    st.divider()

    # Spider Plot section (EBs only — iPSC spider plots not available)
    st.subheader("Lineage Effects (Spider Plot)")

    if perturbations:
        spider_pert = st.selectbox(
            "Select perturbation:",
            options=perturbations,
            key="spider_pert_select",
        )

        spider_path_ebs = get_spider_plot_path(spider_pert, "EBs")
        if spider_path_ebs.exists():
            st.image(str(spider_path_ebs), use_container_width=True)
        else:
            st.info(f"No spider plot for {spider_pert} in EBs")
    else:
        st.warning(f"No perturbations found for {gene}")

    st.divider()

    # Marker Dotplot section
    st.subheader("Marker Gene Expression (Dotplot)")

    if perturbations:
        marker_pert = st.selectbox(
            "Select perturbation for marker expression:",
            options=perturbations,
            key="marker_pert_select",
        )

        st.markdown("**EBs**")
        dotplot_path_ebs = get_marker_dotplot_path(marker_pert, "EBs")
        if dotplot_path_ebs.exists():
            st.image(str(dotplot_path_ebs), use_container_width=True)
        else:
            st.info(f"No marker dotplot for {marker_pert} in EBs")

        st.markdown("**iPSC**")
        dotplot_path_ipsc = get_marker_dotplot_path(marker_pert, "iPSC")
        if dotplot_path_ipsc.exists():
            st.image(str(dotplot_path_ipsc), use_container_width=True)
        else:
            st.info(f"No marker dotplot for {marker_pert} in iPSC")

    st.divider()

    try:
        lineage_df = load_lineage_analysis()
        gene_lineage = lineage_df[lineage_df["gene"] == gene]

        if len(gene_lineage) > 0:
            st.subheader("Lineage Effects")
            fig = plot_lineage_heatmap(gene_lineage, gene)
            st.plotly_chart(fig, use_container_width=True)

            with st.expander("View detailed metrics"):
                display_cols = [
                    "perturbation", "condition", "lineage", "n_cells",
                    "observed_e_score_std", "observed_glass_delta",
                    "q_value", "effect_size_interpretation"
                ]
                available_cols = [c for c in display_cols if c in gene_lineage.columns]
                st.dataframe(gene_lineage[available_cols], use_container_width=True)
        else:
            st.warning(f"No lineage analysis data for {gene}")

    except FileNotFoundError as e:
        st.warning(f"Lineage analysis not available: {e}")

    st.divider()

    # Compositional analysis
    st.subheader("Compositional Changes")

    try:
        comp_df = load_compositional()
        gene_comp = comp_df[comp_df["gene"] == gene]

        if len(gene_comp) > 0:
            fig = plot_compositional_bars(gene_comp, gene)
            st.plotly_chart(fig, use_container_width=True)
        else:
            st.info(f"No significant compositional changes for {gene}")

    except FileNotFoundError as e:
        st.warning(f"Compositional data not available: {e}")

    st.divider()

    # Viability
    st.subheader("Viability / Fitness")

    try:
        via_df = load_viability()

        col1, col2 = st.columns(2)

        with col1:
            fig = plot_viability_rank(via_df, gene, "EBs")
            st.plotly_chart(fig, use_container_width=True)

        with col2:
            fig = plot_viability_rank(via_df, gene, "iPSC")
            st.plotly_chart(fig, use_container_width=True)

    except FileNotFoundError as e:
        st.warning(f"Viability data not available: {e}")

    st.divider()

    # Transcriptome E-distance
    st.subheader("Transcriptome E-distance")

    edist_col1, edist_col2 = st.columns(2)

    with edist_col1:
        try:
            edist_ebs = load_transcriptome_edist("EBs")
            if edist_ebs is not None:
                fig = plot_edist_rank(edist_ebs, gene, "EBs")
                st.plotly_chart(fig, use_container_width=True)
            else:
                st.info("No E-distance data for EBs")
        except Exception as e:
            st.warning(f"E-distance (EBs) not available: {e}")

    with edist_col2:
        try:
            edist_ipsc = load_transcriptome_edist("iPSC")
            if edist_ipsc is not None:
                fig = plot_edist_rank(edist_ipsc, gene, "iPSC")
                st.plotly_chart(fig, use_container_width=True)
            else:
                st.info("No E-distance data for iPSC")
        except Exception as e:
            st.warning(f"E-distance (iPSC) not available: {e}")

    st.divider()

    # Dose Response
    st.subheader("Dose Response")

    dr_col1, dr_col2 = st.columns(2)

    with dr_col1:
        try:
            dr_ebs = load_dose_response("EBs")
            if dr_ebs is not None:
                fig = plot_dose_response_scatter(dr_ebs, gene, "EBs")
                st.plotly_chart(fig, use_container_width=True)
            else:
                st.info("No dose response data for EBs")
        except Exception as e:
            st.warning(f"Dose response (EBs) not available: {e}")

    with dr_col2:
        try:
            dr_ipsc = load_dose_response("iPSC")
            if dr_ipsc is not None:
                fig = plot_dose_response_scatter(dr_ipsc, gene, "iPSC")
                st.plotly_chart(fig, use_container_width=True)
            else:
                st.info("No dose response data for iPSC")
        except Exception as e:
            st.warning(f"Dose response (iPSC) not available: {e}")

    st.divider()

    # TF Similarity
    st.subheader("TF Similarity")

    tf_col1, tf_col2 = st.columns(2)

    for col, cond in [(tf_col1, "EBs"), (tf_col2, "iPSC")]:
        with col:
            st.markdown(f"**{cond}**")
            try:
                tf_df = load_tf_clusters(cond)
                if tf_df is not None:
                    tf_df["gene_name"] = tf_df["perturbation"].str.split("_").str[0]
                    gene_tf = tf_df[tf_df["gene_name"] == gene]
                    if len(gene_tf) > 0:
                        for _, row in gene_tf.iterrows():
                            st.metric(f"{row['perturbation']}", f"Cluster {row['cluster']}")
                    else:
                        st.info(f"No TF cluster data for {gene} in {cond}")

                    # Show overview figure if available
                    tf_fig_path = get_tf_similarity_path(cond)
                    if tf_fig_path.exists():
                        with st.expander("TF Similarity Overview"):
                            st.image(str(tf_fig_path), use_container_width=True)
                else:
                    st.info(f"No TF similarity data for {cond}")
            except Exception as e:
                st.warning(f"TF similarity ({cond}) not available: {e}")

    st.divider()

    # DEG volcano — with condition toggle
    st.subheader("Differential Expression")

    try:
        pert_info = get_perturbations_for_gene(gene)
        perturbations = pert_info["perturbations"]

        if perturbations:
            deg_col1, deg_col2 = st.columns([2, 1])
            with deg_col1:
                selected_pert = st.selectbox(
                    "Select perturbation:",
                    options=perturbations,
                    key="deg_pert_select",
                )
            with deg_col2:
                deg_condition = st.selectbox(
                    "Condition:",
                    options=["EBs", "iPSC"],
                    key="deg_cond_select",
                )

            deg_df = load_deg_table(selected_pert, condition=deg_condition)

            if deg_df is not None:
                fig = plot_deg_volcano(deg_df, f"{selected_pert} ({deg_condition})")
                st.plotly_chart(fig, use_container_width=True)

                st.download_button(
                    "Download DEG table",
                    data=deg_df.to_csv(index=False),
                    file_name=f"{selected_pert}_{deg_condition}_degs.csv",
                    mime="text/csv",
                )

                with st.expander("Top differentially expressed genes"):
                    if "padj" in deg_df.columns:
                        top_degs = deg_df.nsmallest(20, "padj")
                    else:
                        top_degs = deg_df.head(20)
                    st.dataframe(top_degs, use_container_width=True)
            else:
                st.warning(f"No DEG data available for {selected_pert} in {deg_condition}")
        else:
            st.warning("No perturbations found for this gene")

    except Exception as e:
        st.warning(f"DEG data not available: {e}")


# ============ LINEAGE DE TAB ============
with tab_lineage_de:
    st.subheader("Lineage-Specific Differential Expression")

    st.markdown("""
    Differential expression analysis within specific lineages, comparing perturbed vs NTC cells
    that are assigned to the same lineage. Gene IDs are Ensembl (ENSG) identifiers.
    """)

    try:
        pert_info = get_perturbations_for_gene(gene)
        perturbations_lde = pert_info["perturbations"]
    except Exception:
        perturbations_lde = []

    if perturbations_lde:
        lde_col1, lde_col2, lde_col3 = st.columns([2, 1, 1])

        with lde_col1:
            lde_pert = st.selectbox(
                "Select perturbation:",
                options=perturbations_lde,
                key="lde_pert_select",
            )

        with lde_col2:
            lde_cond = st.selectbox(
                "Condition:",
                options=["EBs", "iPSC"],
                key="lde_cond_select",
            )

        # Dynamically list available lineages
        available_lineages = get_available_lineage_de(lde_pert, lde_cond)

        if available_lineages:
            with lde_col3:
                lde_lineage = st.selectbox(
                    "Lineage:",
                    options=available_lineages,
                    key="lde_lineage_select",
                )

            lde_df = load_lineage_de(lde_pert, lde_lineage, lde_cond)

            if lde_df is not None:
                fig = plot_lineage_de_volcano(lde_df, lde_pert, lde_lineage)
                st.plotly_chart(fig, use_container_width=True)

                with st.expander("Top DEGs"):
                    padj_col = "padj" if "padj" in lde_df.columns else "pval"
                    if padj_col in lde_df.columns:
                        top = lde_df.nsmallest(20, padj_col)
                    else:
                        top = lde_df.head(20)
                    st.dataframe(top, use_container_width=True)

                st.download_button(
                    "Download Lineage DE table",
                    data=lde_df.to_csv(index=False),
                    file_name=f"{lde_pert}_{lde_lineage}_{lde_cond}_lineage_de.csv",
                    mime="text/csv",
                )
            else:
                st.warning(f"No lineage DE data for {lde_pert} / {lde_lineage} in {lde_cond}")
        else:
            st.info(f"No lineage DE data available for {lde_pert} in {lde_cond}")
    else:
        st.warning(f"No perturbations found for {gene}")


# ============ PATHWAY ENRICHMENT TAB ============
with tab_pathway:
    st.subheader("Pathway Enrichment Analysis")

    st.markdown("""
    Gene set enrichment analysis of differentially expressed genes per perturbation.
    Shows top enriched pathways/terms ranked by significance.
    """)

    try:
        pert_info = get_perturbations_for_gene(gene)
        perturbations_pe = pert_info["perturbations"]
    except Exception:
        perturbations_pe = []

    if perturbations_pe:
        pe_col1, pe_col2 = st.columns([2, 1])

        with pe_col1:
            pe_pert = st.selectbox(
                "Select perturbation:",
                options=perturbations_pe,
                key="pe_pert_select",
            )

        with pe_col2:
            pe_cond = st.selectbox(
                "Condition:",
                options=["EBs", "iPSC"],
                key="pe_cond_select",
            )

        pe_df = load_pathway_enrichment(pe_cond)

        if pe_df is not None:
            # Filter to selected perturbation
            pe_filtered = pe_df[pe_df["perturbation"] == pe_pert]

            if len(pe_filtered) > 0:
                # Gene set library filter
                if "gene_set_library" in pe_filtered.columns:
                    libraries = ["All"] + sorted(pe_filtered["gene_set_library"].unique().tolist())
                    selected_lib = st.selectbox(
                        "Gene Set Library:",
                        options=libraries,
                        key="pe_lib_select",
                    )
                    if selected_lib != "All":
                        pe_filtered = pe_filtered[pe_filtered["gene_set_library"] == selected_lib]

                fig = plot_pathway_enrichment_bars(pe_filtered, pe_pert)
                st.plotly_chart(fig, use_container_width=True)

                with st.expander("Full Results Table"):
                    st.dataframe(pe_filtered, use_container_width=True)

                st.download_button(
                    "Download Enrichment Results",
                    data=pe_filtered.to_csv(index=False),
                    file_name=f"{pe_pert}_{pe_cond}_enrichment.csv",
                    mime="text/csv",
                )
            else:
                st.info(f"No pathway enrichment data for {pe_pert} in {pe_cond}")
        else:
            st.info(f"No pathway enrichment data available for {pe_cond}")
    else:
        st.warning(f"No perturbations found for {gene}")


# ============ MARKER EXPRESSION TAB ============
with tab_markers:
    st.subheader("Marker / Antibody Expression (TPM)")

    st.markdown("""
    These plots show the expression levels (TPM) of key lineage markers comparing
    **perturbed cells vs NTC (non-targeting control)**. This helps validate that
    lineage composition changes detected in the screen correspond to actual changes
    in marker gene expression.
    """)

    try:
        pert_info = get_perturbations_for_gene(gene)
        perturbations_marker = pert_info["perturbations"]
    except Exception:
        perturbations_marker = []

    if perturbations_marker:
        # Check which perturbations have marker TPM plots available
        available_marker_perts = []
        for p in perturbations_marker:
            for cond in ["EBs", "iPSC"]:
                marker_path = get_marker_tpm_path(p, cond)
                if marker_path.exists():
                    available_marker_perts.append(p)
                    break

        if available_marker_perts:
            marker_pert = st.selectbox(
                "Select perturbation:",
                options=available_marker_perts,
                key="marker_tpm_pert_select",
            )

            # EBs
            marker_path_ebs = get_marker_tpm_path(marker_pert, "EBs")
            if marker_path_ebs.exists():
                st.markdown("**EBs - Marker Expression (TPM)**")
                st.image(str(marker_path_ebs), use_container_width=True)

                st.caption("""
                **How to interpret:** Each bar shows the mean TPM expression of a marker gene.
                Blue = NTC control cells, Orange = perturbed cells.
                Error bars show standard error. Markers are grouped by lineage.
                """)

            # iPSC
            marker_path_ipsc = get_marker_tpm_path(marker_pert, "iPSC")
            if marker_path_ipsc.exists():
                st.divider()
                st.markdown("**iPSC - Marker Expression (TPM)**")
                st.image(str(marker_path_ipsc), use_container_width=True)

            # Antibody validation section
            st.divider()
            st.subheader("Antibody Marker Validation")
            st.markdown("""
            These plots show expression of **antibody-detectable markers** that can be used for
            flow cytometry validation.
            """)

            antibody_path_ebs = get_antibody_validation_path(marker_pert, "EBs")
            antibody_path_ipsc = get_antibody_validation_path(marker_pert, "iPSC")

            has_antibody_ebs = antibody_path_ebs.exists()
            has_antibody_ipsc = antibody_path_ipsc.exists()

            if has_antibody_ebs or has_antibody_ipsc:
                ab_col1, ab_col2 = st.columns(2)

                with ab_col1:
                    st.markdown("**EBs**")
                    if has_antibody_ebs:
                        st.image(str(antibody_path_ebs), use_container_width=True)
                    else:
                        st.info(f"No antibody validation for {marker_pert} in EBs")

                with ab_col2:
                    st.markdown("**iPSC**")
                    if has_antibody_ipsc:
                        st.image(str(antibody_path_ipsc), use_container_width=True)
                    else:
                        st.info(f"No antibody validation for {marker_pert} in iPSC")
            else:
                st.info(f"No antibody validation plots available for {marker_pert}")

        else:
            st.info(f"No marker expression plots available for **{gene}** perturbations.")

            with st.expander("Genes with marker expression data"):
                st.markdown("""
                Marker TPM plots are available for perturbations of these genes:

                ARID2, BCL6, C1orf85, CRTC3, CRX, CSDC2, CSDE1, CTBP2, CTNNB1, EED, EZH2,
                FOXD3, FOXG1, FUS, GATA3, GRHL2, HES5, KDM1A, KDM1B, MLLT1, MLLT10,
                NFKBIE, NKX3-1, NR6A1, NRL, RAX, RLF, SIX2, SOX11, SOX17, SP1, STRAP,
                TADA2B, TBX18, TBX6, TFAP2A, TFAP2B, TFAP2C, TFAP2D, TFAP2E, TGIF1,
                TRIM33, UBTF, USF2, ZBTB12, ZFP2, ZFP90, ZIC2, ZIC3, ZNF177, ZNF200,
                ZNF311, ZNF320, ZNF532, ZNF689, ZNF791
                """)
    else:
        st.warning(f"No perturbations found for {gene}")


# ============ TIMECOURSE TAB ============
with tab_timecourse:
    st.subheader("Expression During Differentiation")

    try:
        tc_df = load_timecourse_expression()
        gene_tc = tc_df[tc_df["gene"] == gene]

        if len(gene_tc) > 0:
            fig = plot_timecourse_expression(tc_df, gene)
            st.plotly_chart(fig, use_container_width=True)

            with st.expander("View timecourse data"):
                st.dataframe(gene_tc, use_container_width=True)
        else:
            st.warning(f"No timecourse data for {gene}")

    except FileNotFoundError as e:
        st.warning(str(e))
        st.info("Run the extraction script to generate timecourse data:")
        st.code("python scripts/extract_timecourse_expression.py")

    st.divider()

    # CellxGene link for timecourse
    st.subheader("Explore in CellxGene")

    try:
        tc_url = st.secrets.get("cellxgene", {}).get("timecourse_url", "")
        if tc_url:
            st.markdown(f"[Open Timecourse in CellxGene]({tc_url})")
        else:
            st.info("CellxGene URL not configured. Add to .streamlit/secrets.toml")
    except Exception:
        st.info("CellxGene URL not configured. Add to .streamlit/secrets.toml")


# ============ CELLXGENE TAB ============
with tab_cellxgene:
    st.subheader("Single-Cell Data Exploration")

    st.markdown("""
    Explore the full single-cell data in CellxGene for detailed analysis:
    """)

    try:
        cellxgene = st.secrets.get("cellxgene", {})

        ipsc_url = cellxgene.get("ipsc_url", "")
        ebs_url = cellxgene.get("ebs_url", "")
        tc_url = cellxgene.get("timecourse_url", "")

        col1, col2, col3 = st.columns(3)

        with col1:
            st.markdown("### iPSC Dataset")
            if ipsc_url:
                st.markdown(f"[Open iPSC in CellxGene]({ipsc_url})")
            else:
                st.warning("URL not configured")

        with col2:
            st.markdown("### EBs Dataset")
            if ebs_url:
                st.markdown(f"[Open EBs in CellxGene]({ebs_url})")
            else:
                st.warning("URL not configured")

        with col3:
            st.markdown("### Timecourse Dataset")
            if tc_url:
                st.markdown(f"[Open Timecourse in CellxGene]({tc_url})")
            else:
                st.warning("URL not configured")

        if not any([ipsc_url, ebs_url, tc_url]):
            st.divider()
            st.info("""
            **To configure CellxGene links:**

            Edit `.streamlit/secrets.toml` and add your CellxGene URLs:

            ```toml
            [cellxgene]
            ipsc_url = "https://your-cellxgene-instance/ipsc"
            ebs_url = "https://your-cellxgene-instance/ebs"
            timecourse_url = "https://your-cellxgene-instance/timecourse"
            ```
            """)

    except Exception as e:
        st.warning(f"Could not load CellxGene configuration: {e}")


# Sidebar
render_sidebar(current_gene=gene)
