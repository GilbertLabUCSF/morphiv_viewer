"""
AI-powered gene summary generation using OpenAI.

Generates detailed, data-driven summaries focused on the MORPHIC TF perturbation screen context.
Includes screen data, lineage effects, DEGs, and mechanistic interpretation.
"""

import os
import streamlit as st
import pandas as pd

try:
    from openai import OpenAI
    OPENAI_AVAILABLE = True
except ImportError:
    OPENAI_AVAILABLE = False

from .data_loader import (
    load_lineage_analysis,
    load_knockdown_efficiency,
    load_viability,
    load_deg_table,
    get_perturbations_for_gene,
)
from .config import LINEAGES

# Note: load_deg_table now accepts an optional condition parameter.
# When called without condition (as below), it searches both EBs and iPSC.


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


def get_gene_screen_data(gene: str) -> dict:
    """
    Extract comprehensive perturbation screen data for a gene.

    Returns dict with knockdown efficiency, viability, lineage effects, and DEGs.
    """
    data = {
        "knockdown": {},
        "viability": {},
        "lineage_effects": {"iPSC": {}, "EBs": {}},
        "top_affected_lineages": [],
        "deg_summary": {},
        "top_degs_up": [],
        "top_degs_down": [],
        "n_cells": {},
        "perturbations": [],
    }

    # Get perturbations for this gene
    try:
        pert_info = get_perturbations_for_gene(gene)
        data["perturbations"] = pert_info.get("perturbations", [])
    except Exception:
        pass

    # Get knockdown efficiency
    try:
        kd_df = load_knockdown_efficiency()
        gene_kd = kd_df[kd_df["target_gene"] == gene]
        for cond in ["iPSC", "EBs"]:
            cond_data = gene_kd[gene_kd["condition"] == cond]
            if len(cond_data) > 0:
                row = cond_data.iloc[0]
                data["knockdown"][cond] = {
                    "pct": round(row["knockdown_pct"], 1),
                    "category": row.get("knockdown_category", ""),
                }
                if "n_cells" in row:
                    data["n_cells"][cond] = int(row["n_cells"])
    except Exception:
        pass

    # Get viability
    try:
        via_df = load_viability()
        gene_via = via_df[via_df["gene"] == gene]
        for cond in ["iPSC", "EBs"]:
            cond_data = gene_via[gene_via["condition"] == cond]
            if len(cond_data) > 0:
                lfc = cond_data.iloc[0]["median_lfc"]
                data["viability"][cond] = {
                    "lfc": round(lfc, 2),
                    "interpretation": get_viability_interpretation(lfc),
                }
    except Exception:
        pass

    # Get lineage effects per condition
    try:
        la_df = load_lineage_analysis()
        gene_la = la_df[la_df["gene"] == gene]

        for cond in ["iPSC", "EBs"]:
            cond_la = gene_la[gene_la["condition"] == cond]
            for lineage in LINEAGES:
                lineage_data = cond_la[cond_la["lineage"] == lineage]
                if len(lineage_data) > 0 and "observed_glass_delta" in lineage_data.columns:
                    deltas = lineage_data["observed_glass_delta"].dropna()
                    if len(deltas) > 0:
                        delta = deltas.iloc[0]
                        data["lineage_effects"][cond][lineage] = round(delta, 2)

        # Get top affected lineages across all conditions
        all_effects = []
        for cond, effects in data["lineage_effects"].items():
            for lineage, delta in effects.items():
                all_effects.append({
                    "condition": cond,
                    "lineage": lineage,
                    "delta": delta,
                    "direction": "enriched" if delta > 0 else "depleted",
                })

        # Sort by absolute delta and take top 5
        all_effects.sort(key=lambda x: abs(x["delta"]), reverse=True)
        data["top_affected_lineages"] = [e for e in all_effects[:5] if abs(e["delta"]) > 0.3]

    except Exception:
        pass

    # Get DEG data
    if data["perturbations"]:
        try:
            # Use first perturbation for DEG data
            deg_df = load_deg_table(data["perturbations"][0])
            if deg_df is not None and len(deg_df) > 0:
                # Count significant DEGs
                if "padj" in deg_df.columns:
                    sig_degs = deg_df[deg_df["padj"] < 0.05]
                    if "log2FoldChange" in deg_df.columns:
                        up_degs = sig_degs[sig_degs["log2FoldChange"] > 0.5]
                        down_degs = sig_degs[sig_degs["log2FoldChange"] < -0.5]
                        data["deg_summary"] = {
                            "total": len(sig_degs),
                            "up": len(up_degs),
                            "down": len(down_degs),
                        }

                        # Get top 5 up and down DEGs
                        gene_col = "gene" if "gene" in deg_df.columns else deg_df.index.name
                        if gene_col and gene_col in deg_df.columns:
                            top_up = up_degs.nlargest(5, "log2FoldChange")
                            top_down = down_degs.nsmallest(5, "log2FoldChange")
                            data["top_degs_up"] = [
                                {"gene": row[gene_col], "lfc": round(row["log2FoldChange"], 2)}
                                for _, row in top_up.iterrows()
                            ]
                            data["top_degs_down"] = [
                                {"gene": row[gene_col], "lfc": round(row["log2FoldChange"], 2)}
                                for _, row in top_down.iterrows()
                            ]
        except Exception:
            pass

    return data


def get_viability_interpretation(lfc: float) -> str:
    """Interpret viability LFC value."""
    if lfc < -1.5:
        return "severely lethal"
    elif lfc < -0.75:
        return "lethal"
    elif lfc < -0.3:
        return "reduced fitness"
    elif lfc > 0.3:
        return "growth advantage"
    else:
        return "neutral"


def format_screen_data_context(gene: str, data: dict) -> str:
    """Format comprehensive screen data as context for the AI prompt."""
    lines = [f"=== MORPHIC CRISPRi Screen Data for {gene} ===\n"]

    # Knockdown efficiency
    if data["knockdown"]:
        lines.append("**Knockdown Efficiency:**")
        for cond, kd in data["knockdown"].items():
            cells = data["n_cells"].get(cond, "?")
            lines.append(f"- {cond}: {kd['pct']}% knockdown ({kd.get('category', '')}) - {cells} cells")
        lines.append("")

    # Viability
    if data["viability"]:
        lines.append("**Viability/Fitness:**")
        for cond, via in data["viability"].items():
            lines.append(f"- {cond}: LFC = {via['lfc']} ({via['interpretation']})")
        lines.append("")

    # Lineage effects
    if data["top_affected_lineages"]:
        lines.append("**Top Lineage Effects (Glass's Delta):**")
        for effect in data["top_affected_lineages"]:
            direction = "↑" if effect["delta"] > 0 else "↓"
            lines.append(f"- {effect['lineage']} ({effect['condition']}): {direction} Δ={effect['delta']} ({effect['direction']})")
        lines.append("")

    # Full lineage breakdown by condition
    for cond in ["EBs", "iPSC"]:
        if data["lineage_effects"].get(cond):
            lines.append(f"**All Lineage Effects in {cond}:**")
            sorted_effects = sorted(data["lineage_effects"][cond].items(), key=lambda x: abs(x[1]), reverse=True)
            for lineage, delta in sorted_effects:
                direction = "enriched" if delta > 0 else "depleted"
                lines.append(f"- {lineage}: Δ={delta} ({direction})")
            lines.append("")

    # DEG summary
    if data["deg_summary"]:
        lines.append(f"**Differential Expression (padj < 0.05, |LFC| > 0.5):**")
        lines.append(f"- Total significant DEGs: {data['deg_summary']['total']}")
        lines.append(f"- Upregulated: {data['deg_summary']['up']}")
        lines.append(f"- Downregulated: {data['deg_summary']['down']}")
        lines.append("")

    # Top DEGs
    if data["top_degs_up"]:
        lines.append("**Top Upregulated Genes:**")
        for deg in data["top_degs_up"]:
            lines.append(f"- {deg['gene']}: LFC = +{deg['lfc']}")
        lines.append("")

    if data["top_degs_down"]:
        lines.append("**Top Downregulated Genes:**")
        for deg in data["top_degs_down"]:
            lines.append(f"- {deg['gene']}: LFC = {deg['lfc']}")
        lines.append("")

    if len(lines) <= 2:
        return ""

    return "\n".join(lines)


GENE_SUMMARY_PROMPT = """You are a scientific expert generating detailed gene summaries for the MORPHIC TF Perturbation Screen Portal.

## About This Screen
This is a CRISPRi screen of transcription factors in human iPSCs and embryoid bodies (EBs), measuring:
- **Lineage composition changes** using Glass's Delta (effect size comparing perturbed vs NTC cells)
- **Viability/fitness** using guide-level log fold changes
- **Transcriptome changes** via differential expression analysis

The lineages tracked are: Amnion, Epiblast, Formative Epiblast, Neural Ectoderm, Non-neural Ectoderm, Trophoblast-Like.

**Effect interpretation:**
- Positive Glass's Delta = lineage ENRICHED after knockdown (more cells in that lineage)
- Negative Glass's Delta = lineage DEPLETED after knockdown (fewer cells in that lineage)
- |Δ| > 0.8 is a strong effect, |Δ| > 0.5 is moderate, |Δ| > 0.3 is small

## Generate a Summary With These Sections:

### Gene Overview
Brief description of {gene}: what type of protein, known functions, expression patterns in development.

### Known Role in Pluripotency & Differentiation
What is established in the literature about this gene's role in:
- Pluripotency maintenance or exit
- Lineage specification (ectoderm, mesoderm, endoderm, trophoblast)
- Chromatin regulation (if applicable)
- Key interacting partners or pathways

### Screen Results Interpretation
Interpret the MORPHIC screen data I provide:
- What do the lineage effects suggest about this gene's function?
- Are the effects consistent with known biology, or surprising?
- What does the DEG pattern (if available) suggest about mechanism?
- Compare iPSC vs EBs effects if both are available.

### Mechanistic Hypothesis
Based on known biology AND the screen results, propose a mechanistic hypothesis:
- What is this gene doing during differentiation?
- Why might knockdown cause the observed lineage changes?
- What downstream targets or pathways might be involved?

### Novel Insights from This Screen
Highlight what this screen reveals that is:
- Confirming known biology
- Potentially novel or unexpected
- Worthy of follow-up validation

### Resources
- [GeneCards](https://www.genecards.org/cgi-bin/carddisp.pl?gene={gene})
- [NCBI Gene](https://www.ncbi.nlm.nih.gov/gene/?term={gene})
- [Human Protein Atlas](https://www.proteinatlas.org/search/{gene})

## Important Guidelines:
- Use bullet points for readability
- Be specific about the screen data - refer to actual numbers
- Connect observations to known biology
- If data suggests something novel, say so explicitly
- Don't cite specific papers, but reference established knowledge
- Aim for ~600-800 words total
"""


@st.cache_data(ttl=86400, show_spinner=False)
def generate_gene_summary(gene_name: str, screen_context: str = "", model: str = "gpt-5.1") -> str:
    """
    Generate an AI-powered gene summary with screen data context.

    Args:
        gene_name: Gene symbol to summarize
        screen_context: Formatted string with screen results for this gene
        model: OpenAI model to use

    Returns:
        Markdown-formatted gene summary
    """
    if not OPENAI_AVAILABLE:
        return None

    client = get_openai_client()
    if client is None:
        return None

    # Build prompt with gene name substitution
    system_prompt = GENE_SUMMARY_PROMPT.replace("{gene}", gene_name)

    # Build user message
    user_message = f"Generate a detailed summary for **{gene_name}**."
    if screen_context:
        user_message += f"\n\n{screen_context}\n\nAnalyze these results in the context of known {gene_name} biology and generate a comprehensive interpretation."

    messages = [
        {"role": "system", "content": system_prompt},
        {"role": "user", "content": user_message}
    ]

    try:
        response = client.chat.completions.create(
            model=model,
            messages=messages,
            temperature=0.5,
            max_tokens=3000,
        )
        return response.choices[0].message.content
    except Exception as e:
        # Fall back to gpt-4o if gpt-5.1 not available
        if "gpt-5" in str(e) or "model" in str(e).lower():
            try:
                response = client.chat.completions.create(
                    model="gpt-4o",
                    messages=messages,
                    temperature=0.5,
                    max_tokens=3000,
                )
                return response.choices[0].message.content
            except Exception as e2:
                return f"Error generating summary: {e2}"
        return f"Error generating summary: {e}"


def render_gene_summary_section(gene: str):
    """
    Render the gene summary section in Streamlit.

    Includes comprehensive screen data context.
    """
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
        """)
        return

    # Get comprehensive screen data for this gene
    screen_data = get_gene_screen_data(gene)
    screen_context = format_screen_data_context(gene, screen_data)

    # Show data availability indicator
    data_available = []
    if screen_data["knockdown"]:
        data_available.append("knockdown")
    if screen_data["viability"]:
        data_available.append("viability")
    if screen_data["top_affected_lineages"]:
        data_available.append("lineage effects")
    if screen_data["deg_summary"]:
        data_available.append(f"{screen_data['deg_summary']['total']} DEGs")

    if data_available:
        st.caption(f"Screen data included: {', '.join(data_available)}")

    # Auto-generate summary (cached for 24 hours)
    with st.spinner(f"Generating AI summary for {gene}..."):
        summary = generate_gene_summary(gene, screen_context)

    if summary:
        # Use expander for cleaner layout
        with st.expander(f"AI-Generated Summary for {gene}", expanded=True):
            st.markdown(summary)
            st.caption("Generated by GPT-5.1 | Cached for 24 hours | Based on MORPHIC screen data + literature")
    else:
        st.error("Failed to generate summary")
