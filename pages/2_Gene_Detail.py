"""
Gene Detail Page - Comprehensive view of a single TF perturbation
"""

import streamlit as st
import pandas as pd

# Add src to path
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent))

from src.data_loader import (
    get_gene_list,
    get_perturbations_for_gene,
    load_lineage_analysis,
    load_knockdown_efficiency,
    load_viability,
    load_compositional,
    load_deg_table,
    load_timecourse_expression,
)
from src.plots import (
    plot_lineage_heatmap,
    plot_knockdown_comparison,
    plot_viability_rank,
    plot_deg_volcano,
    plot_timecourse_expression,
    plot_compositional_bars,
)
from src.config import (
    get_genecards_url,
    get_depmap_url,
    get_ncbi_gene_url,
    get_umap_highlight_path,
    get_marker_dotplot_path,
)
from src.gene_summary import render_gene_summary_section

st.set_page_config(
    page_title="Gene Detail - MORPHIC Portal",
    page_icon=":dna:",
    layout="wide",
)

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
    """Categorize knockdown efficiency: >90% excellent, 70-90% good, 30-70% moderate, <30% poor"""
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
    st.metric("Glass's Δ (max)", f"{glass_delta_max:.2f}" if glass_delta_max is not None else "N/A")

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
tab_overview, tab_perturbation, tab_timecourse, tab_cellxgene = st.tabs([
    "Overview",
    "Perturbation Effects",
    "Timecourse",
    "CellxGene",
])

# ============ OVERVIEW TAB ============
with tab_overview:
    # AI-powered gene summary (general gene info - not perturbation specific)
    render_gene_summary_section(gene)



# ============ PERTURBATION TAB ============
with tab_perturbation:
    # Get perturbation IDs for this gene (used in multiple sections)
    try:
        pert_info = get_perturbations_for_gene(gene)
        perturbations = pert_info["perturbations"]
    except Exception:
        perturbations = []

    # UMAP Highlight section - use pre-generated images from Snakemake pipeline
    st.subheader("Cell Distribution (UMAP)")

    if perturbations:
        # Select perturbation for UMAP
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

    # Marker Dotplot section - use pre-generated images from Snakemake pipeline
    st.subheader("Marker Gene Expression (Dotplot)")

    if perturbations:
        marker_pert = st.selectbox(
            "Select perturbation for marker expression:",
            options=perturbations,
            key="marker_pert_select",
        )

        # EBs dotplot - full width
        st.markdown("**EBs**")
        dotplot_path_ebs = get_marker_dotplot_path(marker_pert, "EBs")
        if dotplot_path_ebs.exists():
            st.image(str(dotplot_path_ebs), use_container_width=True)
        else:
            st.info(f"No marker dotplot for {marker_pert} in EBs")

        # iPSC dotplot - full width
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
            # Lineage heatmap
            st.subheader("Lineage Effects")
            fig = plot_lineage_heatmap(gene_lineage, gene)
            st.plotly_chart(fig, use_container_width=True)

            # Detailed table
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

    # DEG volcano
    st.subheader("Differential Expression")

    # Get perturbation IDs for this gene
    try:
        pert_info = get_perturbations_for_gene(gene)
        perturbations = pert_info["perturbations"]

        if perturbations:
            selected_pert = st.selectbox(
                "Select perturbation:",
                options=perturbations,
                key="deg_pert_select",
            )

            deg_df = load_deg_table(selected_pert)

            if deg_df is not None:
                fig = plot_deg_volcano(deg_df, selected_pert)
                st.plotly_chart(fig, use_container_width=True)

                # Download button
                st.download_button(
                    "Download DEG table",
                    data=deg_df.to_csv(index=False),
                    file_name=f"{selected_pert}_degs.csv",
                    mime="text/csv",
                )

                # Show top DEGs
                with st.expander("Top differentially expressed genes"):
                    if "padj" in deg_df.columns:
                        top_degs = deg_df.nsmallest(20, "padj")
                    else:
                        top_degs = deg_df.head(20)

                    st.dataframe(top_degs, use_container_width=True)
            else:
                st.warning(f"No DEG data available for {selected_pert}")
        else:
            st.warning("No perturbations found for this gene")

    except Exception as e:
        st.warning(f"DEG data not available: {e}")


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
with st.sidebar:
    st.header(f"Gene: {gene}")

    st.caption("Use the tabs above to explore different aspects of this TF perturbation.")

    st.divider()

    st.subheader("External Resources")
    st.markdown(f"[GeneCards]({get_genecards_url(gene)})")
    st.markdown(f"[DepMap]({get_depmap_url(gene)})")
    st.markdown(f"[NCBI Gene]({get_ncbi_gene_url(gene)})")
