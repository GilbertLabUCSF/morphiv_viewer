"""
AI-powered gene summary generation using OpenAI.

Generates GeneCard-style summaries focused on developmental biology context.
"""

import os
import streamlit as st

try:
    from openai import OpenAI
    OPENAI_AVAILABLE = True
except ImportError:
    OPENAI_AVAILABLE = False


def get_openai_client():
    """Get OpenAI client, checking for API key in secrets or environment."""
    api_key = None

    # Try secrets first
    try:
        api_key = st.secrets.get("openai", {}).get("api_key", "")
    except Exception:
        pass

    # Fall back to environment
    if not api_key:
        api_key = os.getenv("OPENAI_API_KEY", "")

    if not api_key:
        return None

    return OpenAI(api_key=api_key)


GENE_SUMMARY_PROMPT = """
You are a knowledgeable scientific assistant who generates GeneCard-style markdown detailed pages.
Be verbose and detailed but don't go overboard.

Your output should include sections such as:
- Overview
- Gene Function in Developmental Biology
- Expression Patterns
- Functional Studies & Perturbation Effects
- Molecular Mechanisms
- Implications for Human Development and Disorders
- Summary of Prior Research
- References & Links

Important requirements:
- NCBI and GeneCards links are important
- Don't cite specific papers
- Pull data from reliable sources like Dependency Map, Human Protein Atlas, Human Cell Atlas and WikiPathways
- Focus on developmental biology context
- Consider roles in pluripotency, lineage commitment, EB vs iPSC differences
- Mention potential effects of knockdown/perturbation

Follow this structure example for POU5F1:

# Gene: POU5F1 (OCT4)

## Overview
**POU5F1**, also known as **OCT4**, is a master transcription factor crucial for maintaining pluripotency...

## Gene Function in Developmental Biology
- **Maintenance of Pluripotency:** ...
- **Cell Fate Determination:** ...

## Expression Patterns
- **Temporal Dynamics:** ...
- **Spatial Distribution:** ...

## Functional Studies & Perturbation Effects
- **Loss-of-Function:** ...
- **Gain-of-Function:** ...

## Molecular Mechanisms
- **DNA-Binding Domains:** ...
- **Transcriptional Network Integration:** ...

## Implications for Human Development and Disorders
- **Developmental Disorders:** ...
- **Cancer and Tumorigenesis:** ...

## Summary of Prior Research
...

## References & Further Reading
1. [NCBI Gene](https://www.ncbi.nlm.nih.gov/gene/...)
2. [GeneCards](https://www.genecards.org/cgi-bin/carddisp.pl?gene=...)
"""


@st.cache_data(ttl=86400, show_spinner=False)
def generate_gene_summary(gene_name: str, model: str = "gpt-4o") -> str:
    """
    Generate an AI-powered gene summary.

    Args:
        gene_name: Gene symbol to summarize
        model: OpenAI model to use (gpt-4o recommended)

    Returns:
        Markdown-formatted gene summary
    """
    if not OPENAI_AVAILABLE:
        return None

    client = get_openai_client()
    if client is None:
        return None

    messages = [
        {"role": "system", "content": GENE_SUMMARY_PROMPT},
        {"role": "user", "content": f"Give me details on {gene_name}"}
    ]

    try:
        response = client.chat.completions.create(
            model=model,
            messages=messages,
            temperature=0.7,
            max_tokens=4096,
            top_p=0.9,
            frequency_penalty=0.3,
            presence_penalty=0.3
        )
        return response.choices[0].message.content
    except Exception as e:
        return f"Error generating summary: {e}"


def render_gene_summary_section(gene: str):
    """
    Render the gene summary section in Streamlit.

    Auto-generates summary on page load with caching.
    """
    st.subheader("AI Gene Summary")

    if not OPENAI_AVAILABLE:
        st.info("OpenAI package not installed. Run: `pip install openai`")
        return

    client = get_openai_client()
    if client is None:
        st.info("""
        **To enable AI gene summaries:**

        Add your OpenAI API key to `.streamlit/secrets.toml`:
        ```toml
        [openai]
        api_key = "sk-..."
        ```

        Or set the `OPENAI_API_KEY` environment variable.
        """)
        return

    # Auto-generate summary (cached for 24 hours)
    with st.spinner(f"Generating AI summary for {gene}..."):
        summary = generate_gene_summary(gene)

    if summary:
        st.markdown(summary)
    else:
        st.error("Failed to generate summary")
