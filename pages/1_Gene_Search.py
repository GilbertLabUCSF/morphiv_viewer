"""
Gene Search Page - Filterable table of all perturbed TFs
"""

import streamlit as st
import pandas as pd

st.set_page_config(
    page_title="Gene Search - MORPHIC Portal",
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

from src.data_loader import load_lineage_analysis, load_knockdown_efficiency, get_perturbations_for_gene
from src.config import CONDITIONS, LINEAGES, EFFECT_SIZES, get_spider_plot_path

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
def sum_unique_cells(group):
    return group.groupby(["perturbation", "condition"])["n_cells"].first().sum()

# Get max absolute Glass's delta per gene
max_glass_delta = (
    filtered_df.groupby("gene")["observed_glass_delta"]
    .apply(lambda x: x.abs().max())
    .reset_index(name="max_glass_delta")
)

# Get lineage effects per gene
def get_lineage_effects(group):
    effects = {}
    for lineage in group["lineage"].unique():
        lineage_data = group[group["lineage"] == lineage]
        deltas = lineage_data["observed_glass_delta"]
        idx_max_abs = deltas.abs().idxmax()
        effects[lineage] = deltas.loc[idx_max_abs]
    return effects

lineage_effects = filtered_df.groupby("gene").apply(get_lineage_effects).reset_index(name="lineage_effects")

summary_df = filtered_df.groupby("gene").agg({
    "perturbation": lambda x: ", ".join(sorted(x.unique())),
    "condition": lambda x: ", ".join(sorted(x.unique())),
}).reset_index()

# Calculate correct cell counts
cell_counts = (
    filtered_df
    .groupby("gene")
    .apply(sum_unique_cells)
    .reset_index(name="n_cells")
)
summary_df = summary_df.merge(cell_counts, on="gene")
summary_df = summary_df.merge(max_glass_delta, on="gene")
summary_df = summary_df.merge(lineage_effects, on="gene")

summary_df.columns = ["Gene", "Perturbation", "Conditions", "Total Cells", "Max |\u0394|", "Lineage Effects"]

# Sort by max glass delta
summary_df = summary_df.sort_values("Max |\u0394|", ascending=False)

# Format lineage effects as colored text
LINEAGE_SHORT = {
    "Amnion": "Amn",
    "Epiblast": "Epi",
    "Formative_Epiblast": "fEpi",
    "Neural_Ectoderm": "NE",
    "Non_neural_Ectoderm": "nnE",
    "Trophoblast_Like": "Troph",
    "Pluripotency": "Pluri",
}

def format_lineage_effects(effects_dict):
    if not effects_dict:
        return ""

    parts = []
    sorted_lineages = sorted(effects_dict.items(), key=lambda x: abs(x[1]), reverse=True)

    for lineage, delta in sorted_lineages:
        short_name = LINEAGE_SHORT.get(lineage, lineage[:4])
        abs_delta = abs(delta)

        if abs_delta < 0.5:
            continue

        if delta > 0:
            parts.append(f"\u25b2{short_name}")
        else:
            parts.append(f"\u25bc{short_name}")

    return " ".join(parts) if parts else "\u2014"

# Display as interactive dataframe with clickable gene links
display_df = summary_df.copy()
display_df["Gene"] = display_df["Gene"].apply(
    lambda g: f"/app/Gene_Detail?gene={g}"
)
display_df["Lineage Effects"] = display_df["Lineage Effects"].apply(format_lineage_effects)

st.dataframe(
    display_df,
    column_config={
        "Gene": st.column_config.LinkColumn(
            "Gene",
            display_text=r"/app/Gene_Detail\?gene=(.+)",
        ),
        "Max |\u0394|": st.column_config.NumberColumn(format="%.2f"),
        "Total Cells": st.column_config.NumberColumn(format="%d"),
    },
    hide_index=True,
    use_container_width=True,
)
st.caption("Click a gene name to view details. Lineage effects: \u25b2 increased, \u25bc decreased (|\u0394| > 0.5)")

# Download button
download_df = summary_df.drop(columns=["Lineage Effects"]).copy()
st.download_button(
    label="Download Filtered Results (CSV)",
    data=download_df.to_csv(index=False),
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

        # Show spider plot preview for selected gene
        try:
            pert_info = get_perturbations_for_gene(selected)
            if pert_info["perturbations"]:
                first_pert = pert_info["perturbations"][0]
                spider_path = get_spider_plot_path(first_pert, "EBs")
                if spider_path.exists():
                    st.image(str(spider_path), width=500, caption=f"Lineage effects for {selected} (EBs)")
        except Exception:
            pass

# Sidebar with quick stats
render_sidebar()

with st.sidebar:
    st.divider()
    st.markdown("**Quick Stats**")
    st.metric("Genes shown", filtered_df["gene"].nunique())
    st.metric("Perturbations", filtered_df["perturbation"].nunique())
