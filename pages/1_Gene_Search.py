"""Search significant lineage phenotypes in the EBs screen."""

from pathlib import Path
import sys
from urllib.parse import quote

import streamlit as st


st.set_page_config(
    page_title="Explore EBs Screen · MORPHIC",
    page_icon="logo_transparent.ico",
    layout="wide",
)

sys.path.insert(0, str(Path(__file__).parent.parent))

from src.components import render_sidebar
from src.config import LINEAGES
from src.data_loader import load_knockdown_efficiency, load_lineage_analysis
from src.data_model import build_gene_results, filter_lineage_hits, lineage_label
from src.styles import inject_css


inject_css()
st.markdown('<p class="eyebrow">EBs PERTURBATION SCREEN</p>', unsafe_allow_html=True)
st.title("Explore lineage phenotypes")
st.markdown(
    "Find transcription-factor knockdowns associated with meaningful shifts in "
    "embryoid-body lineage state. Positive Glass's Δ indicates an increased lineage "
    "score; negative values indicate a decrease."
)

try:
    lineage_df = load_lineage_analysis()
    knockdown_df = load_knockdown_efficiency()
except FileNotFoundError:
    st.error("The EBs screen export is unavailable.")
    st.stop()

st.subheader("Filters")
filter_columns = st.columns([1.15, 1.35, 1, 0.9])
with filter_columns[0]:
    search_term = st.text_input("Gene", placeholder="e.g. GATA3")
with filter_columns[1]:
    lineage_filter = st.multiselect(
        "Lineage",
        options=LINEAGES,
        format_func=lineage_label,
        help="Leave empty to search across every lineage.",
    )
with filter_columns[2]:
    effect_choice = st.selectbox(
        "Minimum effect",
        options=[0.0, 0.3, 0.5, 0.8],
        index=2,
        format_func=lambda value: "Any effect" if value == 0 else f"|Δ| ≥ {value}",
    )
with filter_columns[3]:
    significant_only = st.checkbox(
        "q ≤ 0.05 only",
        value=True,
        help="Use false-discovery-rate-adjusted significance.",
    )

filtered_rows = filter_lineage_hits(
    lineage_df,
    search_term=search_term,
    lineages=lineage_filter,
    min_abs_delta=effect_choice,
    significant_only=significant_only,
)
results = build_gene_results(filtered_rows, knockdown_df)

st.divider()
total_genes = lineage_df["gene"].nunique()
st.subheader(f"{len(results):,} matching genes")
st.caption(
    f"Showing {len(results):,} of {total_genes:,} EBs-screen genes. "
    "Each row reports the strongest lineage among the effects that pass the filters."
)

download_results = results.copy()
display_results = results.copy()
if not display_results.empty:
    display_results["Gene"] = display_results["Gene"].map(
        lambda gene: f"Gene_Detail?gene={quote(str(gene))}"
    )

st.dataframe(
    display_results,
    column_config={
        "Gene": st.column_config.LinkColumn(
            "Gene",
            display_text=r"Gene_Detail\?gene=(.+)",
            help="Open the EBs gene profile",
        ),
        "Max |Δ|": st.column_config.NumberColumn(format="%.2f"),
        "Best q-value": st.column_config.NumberColumn(format="%.2g"),
        "Total cells": st.column_config.NumberColumn(format="%d"),
        "Knockdown": st.column_config.NumberColumn("Median knockdown", format="%.1f%%"),
    },
    column_order=[
        "Gene",
        "Strongest lineage",
        "Effect",
        "Max |Δ|",
        "Best q-value",
        "Knockdown",
        "Perturbations",
        "Total cells",
    ],
    hide_index=True,
    width="stretch",
    height=560,
)

st.download_button(
    "Download these results",
    data=download_results.to_csv(index=False),
    file_name="morphic_ebs_lineage_hits.csv",
    mime="text/csv",
    disabled=download_results.empty,
)

if not results.empty:
    st.markdown("### Open a profile")
    selected_gene = st.selectbox(
        "Selected gene",
        options=results["Gene"].tolist(),
        label_visibility="collapsed",
    )
    st.markdown(
        f'<a class="primary-link" href="Gene_Detail?gene={quote(selected_gene)}">'
        f'Open {selected_gene}&nbsp; →</a>',
        unsafe_allow_html=True,
    )

render_sidebar()

with st.sidebar:
    st.divider()
    st.caption("CURRENT RESULT SET")
    st.metric("Genes", len(results))
    st.metric("Lineage effects", len(filtered_rows))
