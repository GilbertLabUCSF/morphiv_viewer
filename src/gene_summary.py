"""
AI-powered gene summary generation using OpenAI.

Generates data-driven summaries focused on the MORPHIC TF perturbation screen context.
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


def get_openai_client():
    """Get OpenAI client, checking for API key in secrets or environment."""
    api_key = None

    try:
        api_key = st.secrets.get("openai", {}).get("api_key", "")
    except Exception:
        pass

    if not api_key:
        api_key = os.getenv("OPENAI_API_KEY", "")

    if not api_key:
        return None

    return OpenAI(api_key=api_key)


def get_gene_screen_data(gene: str) -> dict:
    """Extract perturbation screen data for a gene."""
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

    try:
        pert_info = get_perturbations_for_gene(gene)
        data["perturbations"] = pert_info.get("perturbations", [])
    except Exception:
        pass

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

        all_effects = []
        for cond, effects in data["lineage_effects"].items():
            for lineage, delta in effects.items():
                all_effects.append({
                    "condition": cond,
                    "lineage": lineage,
                    "delta": delta,
                    "direction": "enriched" if delta > 0 else "depleted",
                })

        all_effects.sort(key=lambda x: abs(x["delta"]), reverse=True)
        data["top_affected_lineages"] = [e for e in all_effects[:5] if abs(e["delta"]) > 0.3]

    except Exception:
        pass

    if data["perturbations"]:
        try:
            deg_df = load_deg_table(data["perturbations"][0])
            if deg_df is not None and len(deg_df) > 0:
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
    """Format screen data as context for the AI prompt."""
    lines = [f"=== MORPHIC CRISPRi Screen Data for {gene} ===\n"]

    if data["knockdown"]:
        lines.append("Knockdown Efficiency:")
        for cond, kd in data["knockdown"].items():
            cells = data["n_cells"].get(cond, "?")
            lines.append(f"- {cond}: {kd['pct']}% ({kd.get('category', '')}) - {cells} cells")
        lines.append("")

    if data["viability"]:
        lines.append("Viability/Fitness:")
        for cond, via in data["viability"].items():
            lines.append(f"- {cond}: LFC = {via['lfc']} ({via['interpretation']})")
        lines.append("")

    if data["top_affected_lineages"]:
        lines.append("Top Lineage Effects (Glass's Delta):")
        for effect in data["top_affected_lineages"]:
            direction = "\u2191" if effect["delta"] > 0 else "\u2193"
            lines.append(f"- {effect['lineage']} ({effect['condition']}): {direction} \u0394={effect['delta']}")
        lines.append("")

    for cond in ["EBs", "iPSC"]:
        if data["lineage_effects"].get(cond):
            lines.append(f"All Lineage Effects in {cond}:")
            sorted_effects = sorted(data["lineage_effects"][cond].items(), key=lambda x: abs(x[1]), reverse=True)
            for lineage, delta in sorted_effects:
                direction = "enriched" if delta > 0 else "depleted"
                lines.append(f"- {lineage}: \u0394={delta} ({direction})")
            lines.append("")

    if data["deg_summary"]:
        lines.append(f"DEGs (padj<0.05, |LFC|>0.5): {data['deg_summary']['total']} total, {data['deg_summary']['up']} up, {data['deg_summary']['down']} down")
        lines.append("")

    if data["top_degs_up"]:
        lines.append("Top Upregulated: " + ", ".join(f"{d['gene']}(+{d['lfc']})" for d in data["top_degs_up"]))
    if data["top_degs_down"]:
        lines.append("Top Downregulated: " + ", ".join(f"{d['gene']}({d['lfc']})" for d in data["top_degs_down"]))

    if len(lines) <= 2:
        return ""

    return "\n".join(lines)


GENE_SUMMARY_PROMPT = """You are a developmental biology expert writing concise gene summaries for the MORPHIC TF Perturbation Screen Portal.

## Screen Context
CRISPRi screen of transcription factors in human iPSCs and embryoid bodies (EBs). Readouts:
- Lineage composition (Glass's Delta): positive = enriched after KD, negative = depleted. |Δ|>0.8 strong, >0.5 moderate, >0.3 small.
- Viability (LFC): negative = fitness cost.
- Transcriptome (DEGs): differentially expressed genes vs NTC.

Lineages: Amnion, Epiblast, Formative Epiblast, Neural Ectoderm, Non-neural Ectoderm, Trophoblast-Like.

## Output Format (4 sections, 300-400 words total)

### What {gene} Does
One paragraph: protein type, known developmental roles, key pathways. Be direct.

### Screen Results
Interpret the data provided. Reference specific numbers. Compare iPSC vs EBs if both available. Bold the most striking finding.

### Mechanistic Hypothesis
2-3 sentences: what is this gene doing during differentiation based on the data? What pathways or targets might explain the lineage shifts?

### Key Takeaways
3-5 bullet points: what confirms known biology, what is novel/unexpected, what deserves follow-up.

## Rules
- Be specific and quantitative — cite the actual Delta values and DEG counts.
- No hedging or filler. State findings directly.
- Do not include resource links (they are shown separately).
- Do not cite specific papers.
"""


@st.cache_data(ttl=86400, show_spinner=False)
def generate_gene_summary(gene_name: str, screen_context: str = "", model: str = "gpt-4.1-mini") -> str:
    """Generate an AI-powered gene summary with screen data context."""
    if not OPENAI_AVAILABLE:
        return None

    client = get_openai_client()
    if client is None:
        return None

    system_prompt = GENE_SUMMARY_PROMPT.replace("{gene}", gene_name)

    user_message = f"Generate a summary for {gene_name}."
    if screen_context:
        user_message += f"\n\n{screen_context}"

    messages = [
        {"role": "system", "content": system_prompt},
        {"role": "user", "content": user_message}
    ]

    try:
        response = client.chat.completions.create(
            model=model,
            messages=messages,
            temperature=0.3,
            max_tokens=1500,
        )
        return response.choices[0].message.content
    except Exception as e:
        # Fall back to gpt-4o-mini
        try:
            response = client.chat.completions.create(
                model="gpt-4o-mini",
                messages=messages,
                temperature=0.3,
                max_tokens=1500,
            )
            return response.choices[0].message.content
        except Exception as e2:
            return f"Error generating summary: {e2}"


def render_gene_summary_section(gene: str):
    """Render the gene summary section in Streamlit."""
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

    screen_data = get_gene_screen_data(gene)
    screen_context = format_screen_data_context(gene, screen_data)

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

    with st.spinner(f"Generating AI summary for {gene}..."):
        summary = generate_gene_summary(gene, screen_context)

    if summary:
        st.markdown(summary)
        st.caption("AI-generated | Cached 24h | Based on screen data + literature")
    else:
        st.error("Failed to generate summary")
