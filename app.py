"""
MORPHIC TF Perturbation Screen Portal - Home Page
"""

import streamlit as st

st.set_page_config(
    page_title="MORPHIC TF Screen Portal",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded",
)

# Custom CSS with Google Fonts
st.markdown("""
<style>
@import url('https://fonts.googleapis.com/css2?family=Inter:wght@300;400;500;600;700&family=JetBrains+Mono:wght@400;500&display=swap');

/* Global font */
html, body, [class*="css"] {
    font-family: 'Inter', -apple-system, BlinkMacSystemFont, sans-serif;
}

/* Headers */
h1, h2, h3, h4, h5, h6 {
    font-family: 'Inter', sans-serif;
    font-weight: 600;
    letter-spacing: -0.02em;
}

/* Main title styling */
.main-title {
    font-size: 2.5rem;
    font-weight: 700;
    letter-spacing: -0.03em;
    margin-bottom: 0.5rem;
    background: linear-gradient(135deg, #60a5fa 0%, #a78bfa 100%);
    -webkit-background-clip: text;
    -webkit-text-fill-color: transparent;
    background-clip: text;
}

/* Subtitle */
.subtitle {
    font-size: 1.1rem;
    color: #94a3b8;
    font-weight: 400;
    margin-bottom: 2rem;
}

/* Metric cards */
.metric-card {
    background: linear-gradient(135deg, rgba(51, 65, 85, 0.5) 0%, rgba(30, 41, 59, 0.5) 100%);
    border: 1px solid rgba(148, 163, 184, 0.1);
    border-radius: 12px;
    padding: 1.25rem;
    text-align: center;
}

.metric-value {
    font-size: 2rem;
    font-weight: 700;
    color: #f1f5f9;
    font-family: 'JetBrains Mono', monospace;
    letter-spacing: -0.02em;
}

.metric-label {
    font-size: 0.85rem;
    color: #94a3b8;
    text-transform: uppercase;
    letter-spacing: 0.05em;
    margin-top: 0.25rem;
}

/* Section headers */
.section-header {
    font-size: 1.25rem;
    font-weight: 600;
    color: #e2e8f0;
    margin: 2rem 0 1rem 0;
    padding-bottom: 0.5rem;
    border-bottom: 1px solid rgba(148, 163, 184, 0.2);
}

/* Info cards */
.info-card {
    background: rgba(51, 65, 85, 0.3);
    border: 1px solid rgba(148, 163, 184, 0.1);
    border-radius: 8px;
    padding: 1rem;
    margin: 0.5rem 0;
}

/* Hide default streamlit branding */
#MainMenu {visibility: hidden;}
footer {visibility: hidden;}
header {visibility: hidden;}

/* Improve dataframe styling */
.stDataFrame {
    font-family: 'Inter', sans-serif;
}

/* Button styling */
.stButton>button {
    font-family: 'Inter', sans-serif;
    font-weight: 500;
}
</style>
""", unsafe_allow_html=True)

# Add src to path
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent))

from src.data_loader import get_gene_list, get_summary_stats

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
left_col, right_col = st.columns([1.5, 1])

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
        ("Lineage Analysis", "Distributional effects on 6 developmental lineages"),
        ("Knockdown Efficiency", "Target gene knockdown validation"),
        ("Viability Scores", "LFC-based cell fitness measurements"),
        ("Differential Expression", "Per-perturbation DEG analysis"),
        ("Timecourse Expression", "Gene expression across differentiation"),
    ]

    for title, desc in data_types:
        st.markdown(f"**{title}** — {desc}")

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
                st.page_link("pages/2_Gene_Detail.py", label="View Gene Details →", icon="🧬")

        else:
            st.warning("No genes found. Check data sources.")

    except Exception as e:
        st.error(f"Error loading gene list: {e}")

    st.markdown("<br>", unsafe_allow_html=True)
    st.markdown('<p class="section-header">Navigation</p>', unsafe_allow_html=True)

    st.page_link("pages/1_Gene_Search.py", label="Gene Search", icon="🔍")
    st.markdown("Filter and explore all perturbed TFs", help=None)

    st.page_link("pages/2_Gene_Detail.py", label="Gene Detail", icon="🧬")
    st.markdown("Detailed view with perturbation effects")

# Footer
st.markdown("---")
st.markdown(
    '<p style="text-align: center; color: #64748b; font-size: 0.85rem;">'
    'MorPhiC Consortium — Molecular Phenotypes of Null Alleles in Cells'
    '</p>',
    unsafe_allow_html=True
)

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

    st.markdown("**Resources**")
    st.link_button("GeneCards", "https://www.genecards.org/", use_container_width=True)
    st.link_button("DepMap", "https://depmap.org/", use_container_width=True)
    st.link_button("NCBI Gene", "https://www.ncbi.nlm.nih.gov/gene/", use_container_width=True)

    st.divider()
    st.caption("v0.1.0")
