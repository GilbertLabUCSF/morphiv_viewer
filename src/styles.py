"""
Shared CSS styles for the MORPHIC Portal.

All pages should call inject_css() to apply consistent styling.
"""

import streamlit as st


SHARED_CSS = """
<style>
@import url('https://fonts.googleapis.com/css2?family=Inter:wght@300;400;500;600;700&display=swap');

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
    font-size: 2.2rem;
    font-weight: 700;
    letter-spacing: -0.03em;
    margin-bottom: 0.25rem;
    background: linear-gradient(135deg, #60a5fa 0%, #a78bfa 100%);
    -webkit-background-clip: text;
    -webkit-text-fill-color: transparent;
    background-clip: text;
}

/* Subtitle */
.subtitle {
    font-size: 1rem;
    color: #94a3b8;
    font-weight: 400;
    margin-bottom: 1.5rem;
}

/* Metric cards */
.metric-card {
    background: linear-gradient(135deg, rgba(51, 65, 85, 0.5) 0%, rgba(30, 41, 59, 0.5) 100%);
    border: 1px solid rgba(148, 163, 184, 0.1);
    border-radius: 12px;
    padding: 1.25rem;
    text-align: center;
    transition: transform 0.15s, box-shadow 0.15s;
}
.metric-card:hover {
    transform: translateY(-2px);
    box-shadow: 0 4px 12px rgba(0, 0, 0, 0.2);
}

.metric-value {
    font-size: 1.75rem;
    font-weight: 700;
    color: #f1f5f9;
    letter-spacing: -0.02em;
}

.metric-label {
    font-size: 0.8rem;
    color: #94a3b8;
    text-transform: uppercase;
    letter-spacing: 0.05em;
    margin-top: 0.25rem;
}

/* Section headers */
.section-header {
    font-size: 1.15rem;
    font-weight: 600;
    color: #e2e8f0;
    margin: 1.5rem 0 0.75rem 0;
    padding-bottom: 0.5rem;
    border-bottom: 1px solid rgba(148, 163, 184, 0.15);
}

/* Figure containers — wrap images/charts for consistent sizing */
.figure-container {
    max-width: 700px;
    margin: 0 auto;
}

/* Hide footer */
footer {visibility: hidden;}

/* Improve dataframe styling */
.stDataFrame {
    font-family: 'Inter', sans-serif;
    font-size: 0.85rem;
}

/* Button styling */
.stButton > button {
    font-family: 'Inter', sans-serif;
    font-weight: 500;
}

/* Better tab styling */
.stTabs [data-baseweb="tab-list"] {
    gap: 4px;
}

.stTabs [data-baseweb="tab"] {
    padding: 8px 16px;
    font-weight: 500;
}

/* Tighter metrics */
[data-testid="stMetric"] {
    background: rgba(51, 65, 85, 0.3);
    border: 1px solid rgba(148, 163, 184, 0.08);
    border-radius: 8px;
    padding: 0.75rem;
}

[data-testid="stMetricLabel"] {
    font-size: 0.75rem !important;
}

[data-testid="stMetricValue"] {
    font-size: 1.25rem !important;
}

/* Constrain image widths in single-column contexts */
[data-testid="stImage"] {
    max-width: 800px;
}

/* Plotly chart max width */
.stPlotlyChart {
    max-width: 900px;
}

/* Better expander styling */
.streamlit-expanderHeader {
    font-size: 0.9rem;
    font-weight: 500;
}

/* Download button less prominent */
.stDownloadButton > button {
    font-size: 0.85rem;
}
</style>
"""


def inject_css():
    """Inject shared CSS styles into the current page."""
    st.markdown(SHARED_CSS, unsafe_allow_html=True)
