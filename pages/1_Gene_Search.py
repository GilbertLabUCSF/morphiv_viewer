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

# Custom CSS with Google Fonts (shared across pages)
st.markdown("""
<style>
@import url('https://fonts.googleapis.com/css2?family=Inter:wght@300;400;500;600;700&family=JetBrains+Mono:wght@400;500&display=swap');

html, body, [class*="css"] {
    font-family: 'Inter', -apple-system, BlinkMacSystemFont, sans-serif;
}

h1, h2, h3, h4, h5, h6 {
    font-family: 'Inter', sans-serif;
    font-weight: 600;
    letter-spacing: -0.02em;
}

.stDataFrame {
    font-family: 'Inter', sans-serif;
}

footer {visibility: hidden;}
</style>
""", unsafe_allow_html=True)

# Add src to path
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent))

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
# Note: n_cells is per perturbation (same value repeated per lineage), so we sum unique values per condition
def sum_unique_cells(group):
    # Get unique n_cells per (perturbation, condition) to avoid counting same cells multiple times
    return group.groupby(["perturbation", "condition"])["n_cells"].first().sum()

# Get max absolute Glass's delta per gene
max_glass_delta = (
    filtered_df.groupby("gene")["observed_glass_delta"]
    .apply(lambda x: x.abs().max())
    .reset_index(name="max_glass_delta")
)

# Get lineage effects per gene (glass delta per lineage, take max across perturbations/conditions)
def get_lineage_effects(group):
    """Get glass delta for each lineage, keeping sign of the value with max absolute magnitude."""
    effects = {}
    for lineage in group["lineage"].unique():
        lineage_data = group[group["lineage"] == lineage]
        # Find the value with the largest absolute magnitude, but keep its sign
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

summary_df.columns = ["Gene", "Perturbation", "Conditions", "Total Cells", "Max |Δ|", "Lineage Effects"]

# Sort by max glass delta
summary_df = summary_df.sort_values("Max |Δ|", ascending=False)

# Format lineage effects as colored text
# Short lineage names for display
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
    """Format lineage effects as text with indicators."""
    if not effects_dict:
        return ""

    parts = []
    # Sort by absolute effect size
    sorted_lineages = sorted(effects_dict.items(), key=lambda x: abs(x[1]), reverse=True)

    for lineage, delta in sorted_lineages:
        short_name = LINEAGE_SHORT.get(lineage, lineage[:4])
        abs_delta = abs(delta)

        # Only show effects with |Δ| > 0.5
        if abs_delta < 0.5:
            continue

        # Simple arrows: ▲ for increased, ▼ for decreased
        if delta > 0:
            parts.append(f"▲{short_name}")
        else:
            parts.append(f"▼{short_name}")

    return " ".join(parts) if parts else "—"

# Display as interactive dataframe with clickable gene links
display_df = summary_df.copy()
display_df["Gene"] = display_df["Gene"].apply(
    lambda g: f"/Gene_Detail?gene={g}"
)
display_df["Lineage Effects"] = display_df["Lineage Effects"].apply(format_lineage_effects)

st.dataframe(
    display_df,
    column_config={
        "Gene": st.column_config.LinkColumn(
            "Gene",
            display_text=r"/Gene_Detail\?gene=(.+)",
        ),
        "Max |Δ|": st.column_config.NumberColumn(format="%.2f"),
        "Total Cells": st.column_config.NumberColumn(format="%d"),
    },
    hide_index=True,
    use_container_width=True,
)
st.caption("Click a gene name to view details. Lineage effects: ▲ increased, ▼ decreased (|Δ| > 0.5)")

# Download button - prepare clean CSV without dict column
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

# Sidebar
with st.sidebar:
    st.markdown("### 🧬 MORPHIC Portal")
    st.caption("TF Perturbation Screen Explorer")

    st.divider()

    st.markdown("**Pages**")
    st.page_link("app.py", label="Home", icon="🏠")
    st.page_link("pages/1_Gene_Search.py", label="Gene Search", icon="🔍")
    st.page_link("pages/2_Gene_Detail.py", label="Gene Detail", icon="🧬")

    st.divider()

    st.markdown("**Quick Stats**")
    st.metric("Genes shown", filtered_df["gene"].nunique())
    st.metric("Perturbations", filtered_df["perturbation"].nunique())

    st.divider()

    st.markdown("**Resources**")
    st.link_button("GeneCards", "https://www.genecards.org/", use_container_width=True)
    st.link_button("DepMap", "https://depmap.org/", use_container_width=True)
    st.link_button("NCBI Gene", "https://www.ncbi.nlm.nih.gov/gene/", use_container_width=True)
