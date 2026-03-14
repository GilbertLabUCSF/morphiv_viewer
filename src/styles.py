"""
Shared CSS styles for the MORPHIC Portal.

All pages should call inject_css() to apply consistent styling.
"""

import streamlit as st


SHARED_CSS = """
<style>
@import url('https://fonts.googleapis.com/css2?family=DM+Serif+Display&family=DM+Sans:wght@300;400;500;600;700&display=swap');

/* Global font */
html, body, [class*="css"] {
    font-family: 'DM Sans', -apple-system, BlinkMacSystemFont, sans-serif;
}

/* Headers */
h1, h2, h3, h4, h5, h6 {
    font-family: 'DM Serif Display', Georgia, serif;
    font-weight: 400;
    letter-spacing: -0.01em;
}

/* Main title styling */
.main-title {
    font-family: 'DM Serif Display', Georgia, serif;
    font-size: 2.2rem;
    font-weight: 400;
    letter-spacing: -0.01em;
    margin-bottom: 0.25rem;
    color: #e8e5e1;
}

/* Subtitle */
.subtitle {
    font-size: 1rem;
    color: #8a8580;
    font-weight: 400;
    margin-bottom: 1.5rem;
}

/* Metric cards */
.metric-card {
    background: #1a1918;
    border: 1px solid rgba(180, 170, 155, 0.1);
    border-radius: 8px;
    padding: 1.25rem;
    text-align: center;
    transition: transform 0.15s, box-shadow 0.15s;
}
.metric-card:hover {
    transform: translateY(-2px);
    box-shadow: 0 4px 12px rgba(0, 0, 0, 0.2);
}

.metric-value {
    font-family: 'DM Serif Display', Georgia, serif;
    font-size: 1.75rem;
    font-weight: 400;
    color: #e8e5e1;
    letter-spacing: -0.02em;
}

.metric-label {
    font-size: 0.8rem;
    color: #8a8580;
    text-transform: uppercase;
    letter-spacing: 0.05em;
    margin-top: 0.25rem;
}

/* Section headers */
.section-header {
    font-family: 'DM Serif Display', Georgia, serif;
    font-size: 1.15rem;
    font-weight: 400;
    color: #e8e5e1;
    margin: 1.5rem 0 0.75rem 0;
    padding-bottom: 0.5rem;
    border-bottom: 1px solid rgba(180, 170, 155, 0.15);
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
    font-family: 'DM Sans', sans-serif;
    font-size: 0.85rem;
}

/* Button styling */
.stButton > button {
    font-family: 'DM Sans', sans-serif;
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
    background: #1a1918;
    border: 1px solid rgba(180, 170, 155, 0.08);
    border-radius: 8px;
    padding: 0.75rem;
}

[data-testid="stMetricLabel"] {
    font-size: 0.75rem !important;
}

[data-testid="stMetricValue"] {
    font-size: 1.25rem !important;
}

/* Subheader spacing */
[data-testid="stSubheader"] {
    margin-top: 0.5rem !important;
    padding-top: 0 !important;
}

/* Tighter dividers */
hr {
    margin: 1rem 0 !important;
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
