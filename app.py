"""MORPHIC EBs TF perturbation browser home page."""

from pathlib import Path
import sys
from urllib.parse import quote

import streamlit as st


st.set_page_config(
    page_title="MORPHIC EBs Perturbation Browser",
    page_icon="logo_transparent.ico",
    layout="wide",
    initial_sidebar_state="auto",
)

sys.path.insert(0, str(Path(__file__).parent))

from src.components import render_sidebar
from src.data_loader import get_data_status, get_gene_list, get_summary_stats
from src.styles import inject_css


inject_css()

st.markdown('<p class="eyebrow">MORPHIC · EMBRYOID BODY SCREEN</p>', unsafe_allow_html=True)
st.markdown('<h1 class="main-title">Explore transcription factor function during early differentiation</h1>', unsafe_allow_html=True)
st.markdown(
    '<p class="subtitle">An EBs-first browser for connecting CRISPRi perturbations '
    'to lineage phenotypes and baseline expression through differentiation.</p>',
    unsafe_allow_html=True,
)

try:
    stats = get_summary_stats()
    if "error" in stats:
        st.error("The primary EBs dataset is unavailable. Please contact the portal maintainers.")
    else:
        metrics = [
            (stats.get("n_genes", 0), "TF genes"),
            (stats.get("n_perturbations", 0), "EB perturbations"),
            (stats.get("n_cells_ebs", 0), "Profiled EB cells"),
            (stats.get("n_timepoints", 0), "Timepoints"),
        ]
        metric_markup = "".join(
            f'<div class="metric-card"><div class="metric-value">{value:,}</div>'
            f'<div class="metric-label">{label}</div></div>'
            for value, label in metrics
        )
        st.markdown(
            f'<div class="metric-grid">{metric_markup}</div>',
            unsafe_allow_html=True,
        )
except Exception:
    st.error("The portal summary could not be loaded.")

st.markdown("<br>", unsafe_allow_html=True)
search_column, context_column = st.columns([1.25, 1], gap="large")

with search_column:
    st.markdown("## Start with a gene")
    st.markdown(
        "Open a profile to see its strongest lineage effects, spider plot, "
        "quality context, and expression trajectory."
    )
    try:
        genes = get_gene_list()
        selected_gene = st.selectbox(
            "Gene",
            options=[""] + genes,
            format_func=lambda value: "Type to search the EBs screen…" if not value else value,
            label_visibility="collapsed",
            key="home_gene_search",
        )
        if selected_gene:
            st.session_state["selected_gene"] = selected_gene
            st.markdown(
                f'<a class="primary-link" href="Gene_Detail?gene={quote(selected_gene)}">'
                f'Open {selected_gene} profile&nbsp; →</a>',
                unsafe_allow_html=True,
            )
        st.page_link(
            "pages/1_Gene_Search.py",
            label="Browse significant lineage effects",
            icon=":material/search:",
        )
    except FileNotFoundError:
        st.error("The EBs gene index is unavailable.")

with context_column:
    st.markdown("## What this browser is for")
    st.markdown(
        """
        - **Find phenotype-shifting TFs** by lineage, direction, effect size, and significance.
        - **Read a perturbation profile** with the spider plot and underlying Glass's Δ values together.
        - **Place the target in developmental context** using unperturbed timecourse expression when present.
        - **Inspect supporting evidence** such as viability, UMAPs, differential expression, and pathways when available.
        """
    )

st.divider()
status = get_data_status()
screen_date = status.get("screen_updated") or "unknown"
timecourse_date = status.get("timecourse_updated") or "not exported"
timecourse_coverage = stats.get("n_screen_genes_timecourse", 0) if "stats" in locals() else 0
screen_gene_count = stats.get("n_genes", 0) if "stats" in locals() else 0
coverage_text = (
    f" It contains expression summaries for **{timecourse_coverage:,} of "
    f"{screen_gene_count:,}** screened targets; absence means the target was not retained "
    "in the processed expression matrix."
    if timecourse_coverage and screen_gene_count
    else ""
)
st.markdown("### Data scope")
st.markdown(
    f"This release focuses on **embryoid bodies (EBs)**. The perturbation screen was "
    f"last exported **{screen_date}**; timecourse summaries were last exported "
    f"**{timecourse_date}**. Timecourse expression is observational baseline context—not "
    f"a perturbation effect.{coverage_text}"
)

st.caption(
    "MorPhiC Consortium · Molecular Phenotypes of Null Alleles in Cells · "
    "Gilbert Lab, UCSF"
)

render_sidebar()
