"""
MORPHIC TF Perturbation Screen Portal - Home Page
"""

import streamlit as st

st.set_page_config(
    page_title="MORPHIC TF Screen Portal",
    page_icon="logo_transparent.ico",
    layout="wide",
    initial_sidebar_state="expanded",
)

# Add src to path
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent))

from src.styles import inject_css
from src.components import render_sidebar
from src.data_loader import get_gene_list, get_summary_stats

inject_css()

# Title with custom styling
st.markdown('<p class="main-title">MORPHIC TF Perturbation Screen</p>', unsafe_allow_html=True)
st.markdown('<p class="subtitle">Exploring transcription factor function in human pluripotent stem cells</p>', unsafe_allow_html=True)

# Summary statistics
try:
    stats = get_summary_stats()

    if "error" in stats:
        st.error(f"Data loading error: {stats['error']}")
    else:
        # Dataset overview metrics
        st.markdown('<p class="section-header">Dataset Overview</p>', unsafe_allow_html=True)

        col1, col2, col3, col4, col5 = st.columns(5)

        metrics = [
            (col1, stats.get('n_genes', 0), "TFs Perturbed"),
            (col2, stats.get('n_perturbations', 0), "Perturbations"),
            (col3, stats.get('n_cells_ipsc', 0), "iPSC Cells"),
            (col4, stats.get('n_cells_ebs', 0), "EBs Cells"),
            (col5, stats.get('n_cells_timecourse', 0), "Timecourse Cells"),
        ]

        for col, value, label in metrics:
            with col:
                st.markdown(f"""
                <div class="metric-card">
                    <div class="metric-value">{value:,}</div>
                    <div class="metric-label">{label}</div>
                </div>
                """, unsafe_allow_html=True)

except Exception as e:
    st.error(f"Error loading summary statistics: {e}")

st.markdown("<br>", unsafe_allow_html=True)

# Two column layout for description and quick search
left_col, right_col = st.columns([3, 2])

with left_col:
    st.markdown('<p class="section-header">About</p>', unsafe_allow_html=True)

    st.markdown("""
    The **MorPhiC Consortium** (Molecular Phenotypes of Null Alleles in Cells) is an NIH
    initiative to catalogue molecular and cellular phenotypes from inactivating every human gene.

    This portal provides access to the **TF Perturbation Screen**, with CRISPRi perturbations
    of transcription factors in iPSCs and embryoid bodies (EBs).
    """)

    st.markdown('<p class="section-header">Available Data</p>', unsafe_allow_html=True)

    data_types = [
        ("Lineage Analysis", "Effects on 6 developmental lineages"),
        ("Knockdown Efficiency", "Target gene knockdown validation"),
        ("Viability", "Cell fitness measurements"),
        ("Differential Expression", "DEG analysis (EBs & iPSC)"),
        ("Timecourse", "Expression across differentiation"),
        ("Dose Response", "Knockdown vs lineage effect"),
        ("Lineage DE", "Within-lineage DE (perturbed vs NTC)"),
        ("Pathway Enrichment", "Gene set enrichment of DEGs"),
        ("TF Similarity", "TF clustering by transcriptomic signatures"),
        ("Transcriptome E-distance", "Global transcriptome shift"),
    ]

    # Two-column data type list
    dt_left, dt_right = st.columns(2)
    for i, (title, desc) in enumerate(data_types):
        with dt_left if i % 2 == 0 else dt_right:
            st.markdown(f"**{title}** \u2014 {desc}")

with right_col:
    st.markdown('<p class="section-header">Quick Gene Search</p>', unsafe_allow_html=True)

    try:
        genes = get_gene_list()

        if genes:
            selected_gene = st.selectbox(
                "Search for a gene:",
                options=[""] + genes,
                format_func=lambda x: "Type to search..." if x == "" else x,
                key="home_gene_search",
                label_visibility="collapsed",
            )

            if selected_gene:
                st.session_state["selected_gene"] = selected_gene
                st.success(f"**{selected_gene}** selected")
                st.page_link("pages/2_Gene_Detail.py", label="View Gene Details \u2192", icon="\U0001f9ec")

        else:
            st.warning("No genes found. Check data sources.")

    except Exception as e:
        st.error(f"Error loading gene list: {e}")

    st.markdown("<br>", unsafe_allow_html=True)
    st.markdown('<p class="section-header">Navigation</p>', unsafe_allow_html=True)

    st.page_link("pages/1_Gene_Search.py", label="Gene Search", icon="\U0001f50d")
    st.markdown("Filter and explore all perturbed TFs", help=None)

    st.page_link("pages/2_Gene_Detail.py", label="Gene Detail", icon="\U0001f9ec")
    st.markdown("Detailed view with perturbation effects")

# Footer
st.markdown("---")
st.markdown(
    '<p style="text-align: center; color: #64748b; font-size: 0.85rem;">'
    'MorPhiC Consortium \u2014 Molecular Phenotypes of Null Alleles in Cells'
    '</p>',
    unsafe_allow_html=True
)

# Sidebar
render_sidebar()
