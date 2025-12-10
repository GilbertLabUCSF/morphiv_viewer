"""
MORPHIC TF Perturbation Screen Portal - Landing Page
"""

import streamlit as st

st.set_page_config(
    page_title="MORPHIC TF Screen Portal",
    page_icon=":dna:",
    layout="wide",
    initial_sidebar_state="expanded",
)

# Add src to path
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent))

from src.data_loader import get_gene_list, get_summary_stats

# Title
st.title("MORPHIC Transcription Factor Perturbation Screen")

st.markdown("""
The **MorPhiC Consortium** (Molecular Phenotypes of Null Alleles in Cells) is a National Institutes of Health (NIH)
initiative to create a comprehensive catalogue of the molecular and cellular phenotypes resulting from the
inactivation of every human gene.

This portal provides access to the **TF Perturbation Screen** data, enabling exploration of transcription factors
perturbed via CRISPRi in human pluripotent stem cells (iPSCs) and embryoid bodies (EBs).
""")

# Summary statistics
st.header("Dataset Overview")

try:
    stats = get_summary_stats()

    if "error" in stats:
        st.error(f"Data loading error: {stats['error']}")
    else:
        col1, col2, col3, col4, col5 = st.columns(5)

        with col1:
            st.metric("TFs Perturbed", f"{stats.get('n_genes', 'N/A'):,}")

        with col2:
            st.metric("Perturbations", f"{stats.get('n_perturbations', 'N/A'):,}")

        with col3:
            st.metric("iPSC Cells", f"{stats.get('n_cells_ipsc', 0):,}")

        with col4:
            st.metric("EBs Cells", f"{stats.get('n_cells_ebs', 0):,}")

        with col5:
            st.metric("Timecourse Cells", f"{stats.get('n_cells_timecourse', 0):,}")

except Exception as e:
    st.error(f"Error loading summary statistics: {e}")

st.divider()

# Quick gene search
st.header("Quick Gene Search")

try:
    genes = get_gene_list()

    if genes:
        selected_gene = st.selectbox(
            "Search for a gene:",
            options=[""] + genes,
            format_func=lambda x: "Type to search..." if x == "" else x,
            key="home_gene_search",
        )

        if selected_gene:
            st.info(f"Navigate to **Gene Detail** page in the sidebar to view {selected_gene}")
            # Store selected gene in session state for the detail page
            st.session_state["selected_gene"] = selected_gene

    else:
        st.warning("No genes found. Check that data sources are properly linked.")

except Exception as e:
    st.error(f"Error loading gene list: {e}")

st.divider()

# About section
st.header("About This Portal")

st.markdown("""
### Data Available

- **Lineage Analysis**: Distributional effects on 6 developmental lineages
  - Amnion, Epiblast, Formative Epiblast, Neural Ectoderm, Non-neural Ectoderm, Trophoblast-Like

- **Knockdown Efficiency**: Target gene knockdown validation

- **Viability/Fitness**: LFC-based cell fitness scores

- **Differential Expression**: Per-perturbation DEG analysis

- **Timecourse Expression**: Gene expression across differentiation

### Navigation

Use the sidebar to navigate between pages:

1. **Gene Search**: Filter and explore all perturbed TFs
2. **Gene Detail**: Detailed view of individual genes with perturbation effects,
   timecourse expression, and links to CellxGene

### External Resources

- [GeneCards](https://www.genecards.org/) - Gene annotations
- CellxGene - Single-cell data exploration (links on gene pages)

---

*MorPhiC Consortium - Molecular Phenotypes of Null Alleles in Cells*
""")

# Sidebar info
with st.sidebar:
    st.header("MORPHIC Portal")
    st.caption("TF Perturbation Screen Explorer")

    st.divider()

    st.markdown("**Pages**")
    st.page_link("pages/1_Gene_Search.py", label="Gene Search", icon="🔍")
    st.page_link("pages/2_Gene_Detail.py", label="Gene Detail", icon="🧬")

    st.divider()

    st.caption("v0.1.0 - MVP")
