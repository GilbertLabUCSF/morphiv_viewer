"""Pure data transformations used by the EBs portal UI."""

from __future__ import annotations

import pandas as pd


def lineage_label(value: str) -> str:
    """Format pipeline lineage identifiers for display."""
    labels = {
        "Amnion": "Amnion",
        "Epiblast": "Epiblast",
        "Formative_Epiblast": "Formative epiblast",
        "Neural_Ectoderm": "Neural ectoderm",
        "Non_neural_Ectoderm": "Non-neural ectoderm",
        "Trophoblast_Like": "Trophoblast-like",
    }
    return labels.get(str(value), str(value).replace("_", " "))


def timepoint_label(value: object) -> str:
    """Format developmental timepoints without over-emphasizing cell provenance."""
    text = str(value)
    return "Baseline" if text == "iPSC" else f"Day {text}"


def filter_lineage_hits(
    lineage_df: pd.DataFrame,
    search_term: str = "",
    lineages: list[str] | None = None,
    min_abs_delta: float = 0.5,
    significant_only: bool = True,
) -> pd.DataFrame:
    """Filter lineage rows using literal gene matching and real effect criteria."""
    filtered = lineage_df.copy()

    if search_term:
        filtered = filtered[
            filtered["gene"].str.contains(
                search_term.strip(), case=False, na=False, regex=False
            )
        ]

    if lineages:
        filtered = filtered[filtered["lineage"].isin(lineages)]

    if min_abs_delta > 0:
        filtered = filtered[
            filtered["observed_glass_delta"].abs() >= min_abs_delta
        ]

    if significant_only and "q_value" in filtered.columns:
        filtered = filtered[filtered["q_value"] <= 0.05]

    return filtered


def build_gene_results(
    filtered_lineage_df: pd.DataFrame,
    knockdown_df: pd.DataFrame | None = None,
) -> pd.DataFrame:
    """Build one concise result row per gene from filtered lineage effects."""
    columns = [
        "Gene",
        "Strongest lineage",
        "Effect",
        "Max |Δ|",
        "Best q-value",
        "Perturbations",
        "Total cells",
        "Knockdown",
    ]
    if filtered_lineage_df.empty:
        return pd.DataFrame(columns=columns)

    df = filtered_lineage_df.copy()
    df["abs_delta"] = df["observed_glass_delta"].abs()
    strongest = df.loc[df.groupby("gene")["abs_delta"].idxmax()].copy()
    strongest = strongest.set_index("gene")

    perturbations = (
        df.groupby("gene")["perturbation"]
        .agg(lambda values: ", ".join(sorted(values.unique())))
    )
    cell_counts = (
        df.drop_duplicates(["gene", "perturbation"])
        .groupby("gene")["n_cells"]
        .sum()
    )
    best_q = df.groupby("gene")["q_value"].min()

    result = pd.DataFrame(index=strongest.index)
    result["Strongest lineage"] = strongest["lineage"].map(lineage_label)
    result["Effect"] = strongest["observed_glass_delta"].map(
        lambda value: "Increased" if value > 0 else "Decreased"
    )
    result["Max |Δ|"] = strongest["abs_delta"]
    result["Best q-value"] = best_q
    result["Perturbations"] = perturbations
    result["Total cells"] = cell_counts

    if knockdown_df is not None and not knockdown_df.empty:
        knockdown = knockdown_df.groupby("target_gene")["knockdown_pct"].median()
        result["Knockdown"] = result.index.to_series().map(knockdown)
    else:
        result["Knockdown"] = pd.NA

    result.index.name = "Gene"
    result = result.reset_index()
    return result[columns].sort_values("Max |Δ|", ascending=False)


def perturbation_profile(lineage_df: pd.DataFrame, perturbation: str) -> dict:
    """Return headline metrics for one perturbation."""
    rows = lineage_df[lineage_df["perturbation"] == perturbation].copy()
    if rows.empty:
        return {}

    rows["abs_delta"] = rows["observed_glass_delta"].abs()
    strongest = rows.loc[rows["abs_delta"].idxmax()]
    significant = rows[rows["q_value"] <= 0.05] if "q_value" in rows else rows.iloc[0:0]
    return {
        "n_cells": int(strongest["n_cells"]),
        "strongest_lineage": lineage_label(strongest["lineage"]),
        "strongest_delta": float(strongest["observed_glass_delta"]),
        "strongest_q": float(strongest["q_value"]),
        "n_significant_lineages": int(significant["lineage"].nunique()),
        "significant_effects": [
            {
                "lineage": lineage_label(row["lineage"]),
                "delta": float(row["observed_glass_delta"]),
                "q_value": float(row["q_value"]),
            }
            for _, row in significant.sort_values(
                "abs_delta", ascending=False
            ).iterrows()
        ],
    }
