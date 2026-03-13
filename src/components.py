"""
Shared UI components for the MORPHIC Portal.

Provides reusable sidebar and other common elements.
"""

import streamlit as st


def render_sidebar(current_gene: str = None):
    """Render the shared sidebar navigation.

    Args:
        current_gene: If set, shows gene-specific links in the sidebar.
    """
    with st.sidebar:
        st.markdown("### MORPHIC Portal")
        st.caption("TF Perturbation Screen Explorer")

        st.divider()

        st.markdown("**Pages**")
        st.page_link("app.py", label="Home", icon="\U0001f3e0")
        st.page_link("pages/1_Gene_Search.py", label="Gene Search", icon="\U0001f50d")
        st.page_link("pages/2_Gene_Detail.py", label="Gene Detail", icon="\U0001f9ec")

        st.divider()

        if current_gene:
            from .config import get_genecards_url, get_depmap_url, get_ncbi_gene_url

            st.markdown(f"**Current Gene: {current_gene}**")
            st.caption("Use tabs to explore perturbation data")

            st.divider()

            st.markdown(f"**Resources for {current_gene}**")
            st.link_button(f"GeneCards: {current_gene}", get_genecards_url(current_gene), use_container_width=True)
            st.link_button(f"DepMap: {current_gene}", get_depmap_url(current_gene), use_container_width=True)
            st.link_button(f"NCBI Gene: {current_gene}", get_ncbi_gene_url(current_gene), use_container_width=True)
        else:
            st.markdown("**Resources**")
            st.link_button("GeneCards", "https://www.genecards.org/", use_container_width=True)
            st.link_button("DepMap", "https://depmap.org/", use_container_width=True)
            st.link_button("NCBI Gene", "https://www.ncbi.nlm.nih.gov/gene/", use_container_width=True)

            st.divider()
            st.caption("v0.2.0")
