"""
Gene Search Page - Filterable table of all perturbed TFs
"""

import streamlit as st
import pandas as pd

# Add src to path
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent))

from src.data_loader import load_lineage_analysis, load_knockdown_efficiency
from src.config import CONDITIONS, LINEAGES, EFFECT_SIZES

st.set_page_config(
    page_title="Gene Search - MORPHIC Portal",
    page_icon=":mag:",
    layout="wide",
)

st.title("Gene Search")
st.markdown("Search and filter perturbed transcription factors.")

# Load data
try:
    lineage_df = load_lineage_analysis()
    knockdown_df = load_knockdown_efficiency()
except FileNotFoundError as e:
    st.error(f"Data not available: {e}")
    st.stop()

# Filters in columns
st.subheader("Filters")

col1, col2, col3, col4 = st.columns(4)

with col1:
    search_term = st.text_input("Gene name:", placeholder="e.g., SOX2")

with col2:
    condition_filter = st.multiselect(
        "Condition:",
        options=CONDITIONS,
        default=CONDITIONS,
    )

with col3:
    lineage_filter = st.multiselect(
        "Lineage:",
        options=LINEAGES,
        default=[],
        help="Filter to genes affecting these lineages",
    )

with col4:
    effect_filter = st.selectbox(
        "Effect size:",
        options=["All"] + EFFECT_SIZES,
    )

# Apply filters
filtered_df = lineage_df.copy()

# Gene name search
if search_term:
    filtered_df = filtered_df[
        filtered_df["gene"].str.contains(search_term, case=False, na=False)
    ]

# Condition filter
if condition_filter:
    filtered_df = filtered_df[filtered_df["condition"].isin(condition_filter)]

# Lineage filter
if lineage_filter:
    filtered_df = filtered_df[filtered_df["lineage"].isin(lineage_filter)]

# Effect size filter
if effect_filter != "All":
    filtered_df = filtered_df[filtered_df["effect_size_interpretation"] == effect_filter]

st.divider()

# Summary of filtered results
st.subheader(f"Results: {filtered_df['gene'].nunique()} genes")

# Create summary table (one row per gene)
# Note: n_cells is per perturbation (same value repeated per lineage), so we sum unique values per condition
def sum_unique_cells(group):
    # Get unique n_cells per (perturbation, condition) to avoid counting same cells multiple times
    return group.groupby(["perturbation", "condition"])["n_cells"].first().sum()

summary_df = filtered_df.groupby("gene").agg({
    "perturbation": lambda x: ", ".join(sorted(x.unique())),
    "condition": lambda x: ", ".join(sorted(x.unique())),
    "max_abs_z_score": "max",
    "lineage": lambda x: ", ".join(sorted(x.unique())),
}).reset_index()

# Calculate correct cell counts
cell_counts = (
    filtered_df
    .groupby("gene")
    .apply(sum_unique_cells)
    .reset_index(name="n_cells")
)
summary_df = summary_df.merge(cell_counts, on="gene")

summary_df.columns = ["Gene", "Perturbation", "Conditions", "Max Z-Score", "Lineages Tested", "Total Cells"]

# Sort by max z-score
summary_df = summary_df.sort_values("Max Z-Score", ascending=False)

# Display as interactive dataframe with clickable gene links
display_df = summary_df.rename(columns={"Max Z-Score": "Max |Z|"}).copy()
display_df["Gene"] = display_df["Gene"].apply(
    lambda g: f"/Gene_Detail?gene={g}"
)
st.dataframe(
    display_df,
    column_config={
        "Gene": st.column_config.LinkColumn(
            "Gene",
            display_text=r"/Gene_Detail\?gene=(.+)",
        ),
        "Max |Z|": st.column_config.NumberColumn(format="%.2f"),
        "Total Cells": st.column_config.NumberColumn(format="%d"),
    },
    hide_index=True,
    use_container_width=True,
)
st.caption("Click a gene name to view details.")

# Download button
st.download_button(
    label="Download Filtered Results (CSV)",
    data=summary_df.to_csv(index=False),
    file_name="morphic_gene_search_results.csv",
    mime="text/csv",
)

st.divider()

# Click to view gene detail
st.subheader("View Gene Details")

if len(summary_df) > 0:
    selected = st.selectbox(
        "Select a gene to view details:",
        options=[""] + summary_df["Gene"].tolist(),
        format_func=lambda x: "Choose a gene..." if x == "" else x,
    )

    if selected:
        st.session_state["selected_gene"] = selected
        st.info(f"Gene **{selected}** selected. Go to **Gene Detail** page in the sidebar.")

# Sidebar
with st.sidebar:
    st.header("Quick Stats")

    st.metric("Genes shown", filtered_df["gene"].nunique())
    st.metric("Total perturbations", filtered_df["perturbation"].nunique())

    st.divider()

    st.caption("Tip: Use filters to narrow down results, then click a gene to view details.")
