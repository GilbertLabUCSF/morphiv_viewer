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

/* Use the curated navigation below, not Streamlit's duplicate page list. */
[data-testid="stSidebarNav"] { display: none; }

.eyebrow {
    color: #c8a57b;
    font-size: 0.72rem;
    font-weight: 700;
    letter-spacing: 0.12em;
    margin: 0 0 0.65rem 0;
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
    font-size: clamp(2.2rem, 4vw, 3.5rem);
    font-weight: 400;
    letter-spacing: -0.01em;
    line-height: 1.08;
    max-width: 980px;
    margin-bottom: 0.65rem;
    color: #e8e5e1;
}

/* Subtitle */
.subtitle {
    font-size: 1.08rem;
    color: #8a8580;
    font-weight: 400;
    line-height: 1.65;
    max-width: 820px;
    margin-bottom: 1.75rem;
}

.primary-link {
    display: inline-flex;
    align-items: center;
    margin: 0.65rem 0 1rem;
    padding: 0.72rem 1rem;
    border-radius: 7px;
    background: #c8a57b;
    color: #0c0b0a !important;
    font-weight: 700;
    text-decoration: none !important;
}
.primary-link:hover { background: #d4b48e; }

.evidence-callout {
    margin: 0.45rem 0 0.25rem;
    padding: 0.72rem 0.9rem;
    border-left: 3px solid #c8a57b;
    border-radius: 4px;
    background: rgba(200, 165, 123, 0.08);
    color: #d7d2cc;
    line-height: 1.45;
}

.profile-heading {
    display: flex;
    align-items: baseline;
    gap: 1rem;
    margin: 0.55rem 0 0.7rem;
}
.profile-heading span {
    font-family: 'DM Serif Display', Georgia, serif;
    font-size: 2rem;
    color: #e8e5e1;
}
.profile-heading small { color: #8a8580; }
.profile-heading small::before { content: "·"; margin-right: 1rem; }

/* Metric cards */
.metric-grid {
    display: grid;
    grid-template-columns: repeat(4, minmax(0, 1fr));
    gap: 1rem;
}

.metric-card {
    background: #1a1918;
    border: 1px solid rgba(180, 170, 155, 0.1);
    border-radius: 8px;
    padding: 1.25rem;
    text-align: center;
    transition: transform 0.15s, box-shadow 0.15s;
}

.profile-metric-grid {
    display: grid;
    grid-template-columns: repeat(5, minmax(0, 1fr));
    gap: 0.8rem;
}
.profile-metric-grid .metric-card {
    min-height: 104px;
    padding: 0.9rem 1rem;
    text-align: left;
}
.profile-metric-label { color: #c9c4bd; font-size: 0.83rem; }
.profile-metric-value { color: #f1eee9; font-size: 1.38rem; margin-top: 0.25rem; }
.profile-metric-note { color: #8a8580; font-size: 0.76rem; margin-top: 0.15rem; }
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

/* Keep secondary text readable on dark surfaces. */
.stCaption, [data-testid="stCaptionContainer"] { color: #8a8580 !important; }

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

@media (max-width: 700px) {
    .main-title { font-size: 2.15rem; }
    .metric-grid, .profile-metric-grid {
        grid-template-columns: repeat(2, minmax(0, 1fr));
        gap: 0.55rem;
    }
    .metric-card { padding: 0.9rem 0.5rem; }
    .metric-value { font-size: 1.4rem; }
    .profile-heading { display: block; }
    .profile-heading small { display: block; margin-top: 0.2rem; }
    .profile-heading small::before { content: ""; margin: 0; }
    .stTabs [data-baseweb="tab"] { padding: 7px 9px; }
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
