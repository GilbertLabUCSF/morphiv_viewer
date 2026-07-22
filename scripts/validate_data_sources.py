#!/usr/bin/env python3
"""Validate the EBs-focused portal export before startup or deployment."""

from __future__ import annotations

import os
from pathlib import Path
import sys

import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parent.parent
DATA_ROOT = Path(os.environ.get("MORPHIC_DATA_PATH", PROJECT_ROOT / "data"))
EXTRACTED_ROOT = PROJECT_ROOT / "data_extracted"
RESULTS = DATA_ROOT / "tf_perturbseq" / "results" / "EBs"
FIGURES = DATA_ROOT / "tf_perturbseq" / "figures" / "EBs"

CORE_TABLES = {
    "lineage": (
        RESULTS / "lineage_analysis" / "EBs_lineage_analysis_merged.csv",
        {"perturbation", "lineage", "n_cells", "observed_glass_delta", "q_value"},
    ),
    "knockdown": (
        RESULTS / "knockdown_efficiency" / "knockdown_efficiency_all_genes.csv",
        {"perturbation", "target_gene", "knockdown_pct", "n_cells"},
    ),
    "viability": (
        RESULTS / "viability" / "viability_scores_gene_level.csv",
        {"gene", "median_lfc", "stouffer_q"},
    ),
    "probability shifts": (
        RESULTS / "knn_probability_shifts" / "EBs_knn_prob_shifts.csv",
        {
            "perturbation",
            "lineage",
            "mean_prob_pert",
            "mean_prob_ntc",
            "log2fc",
            "ci_low",
            "ci_high",
            "q_value_boot",
            "n_cells",
        },
    ),
}

OPTIONAL_ASSETS = {
    "Compositional hits": RESULTS / "compositional" / "EBs_significant_hits.csv",
    "Transcriptome E-distance": RESULTS / "transcriptome_edist" / "EBs_tedist_merged.csv",
    "Dose response": RESULTS / "dose_response" / "EBs_dose_response.csv",
    "TF clusters": RESULTS / "tf_similarity" / "EBs_tf_clusters.csv",
    "Pathway enrichment lookup": EXTRACTED_ROOT / "pathway_enrichment_ebs.parquet",
    "DEG tables": RESULTS / "differential_expression" / "deg_tables",
    "Lineage DE": RESULTS / "lineage_de",
    "Spider plots": FIGURES / "lineage_analysis" / "spider",
    "UMAP highlights": FIGURES / "umap" / "highlights",
}

TIMECOURSE_SOURCE = (
    DATA_ROOT
    / "single_cell_timecourse"
    / "results"
    / "processed"
    / "timecourse_merged.h5ad"
)
TIMECOURSE_EXPORT = EXTRACTED_ROOT / "timecourse_expression.parquet"


def validate_table(label: str, path: Path, required_columns: set[str]) -> list[str]:
    issues = []
    if not path.exists():
        return [f"{label}: missing {path}"]
    try:
        sample = pd.read_csv(path, nrows=10)
    except Exception as exc:
        return [f"{label}: unreadable ({exc})"]
    missing_columns = required_columns - set(sample.columns)
    if missing_columns:
        issues.append(f"{label}: missing columns {sorted(missing_columns)}")
    if path.stat().st_size == 0:
        issues.append(f"{label}: file is empty")
    return issues


def main() -> int:
    print("MORPHIC EBs portal data validation")
    print("=" * 40)
    issues: list[str] = []

    if not DATA_ROOT.exists():
        issues.append(f"Data root is unavailable: {DATA_ROOT}")
    else:
        print(f"✓ Data root: {DATA_ROOT.resolve()}")

    for label, (path, columns) in CORE_TABLES.items():
        table_issues = validate_table(label, path, columns)
        issues.extend(table_issues)
        if not table_issues:
            print(f"✓ Core table: {label}")

    if not TIMECOURSE_EXPORT.exists():
        issues.append(
            "Timecourse export is missing; run scripts/extract_timecourse_expression.py"
        )
    else:
        try:
            timecourse = pd.read_parquet(TIMECOURSE_EXPORT)
            required = {
                "gene",
                "day",
                "mean_expression",
                "pct_detected",
                "n_cells",
                "expression_scale",
            }
            missing = required - set(timecourse.columns)
            if missing:
                issues.append(f"Timecourse export is missing columns {sorted(missing)}")
            elif timecourse.empty:
                issues.append("Timecourse export is empty")
            else:
                target_path = CORE_TABLES["knockdown"][0]
                coverage_text = ""
                if target_path.exists():
                    targets = pd.read_csv(target_path, usecols=["target_gene"])
                    target_genes = set(targets["target_gene"].dropna().astype(str))
                    timecourse_genes = set(timecourse["gene"].dropna().astype(str))
                    covered_targets = len(target_genes & timecourse_genes)
                    coverage_text = (
                        f", {covered_targets:,}/{len(target_genes):,} "
                        "screen targets represented"
                    )
                print(
                    f"✓ Timecourse export: {timecourse['gene'].nunique():,} genes, "
                    f"{timecourse['day'].nunique()} timepoints{coverage_text}"
                )
        except Exception as exc:
            issues.append(f"Timecourse export is unreadable ({exc})")

    if TIMECOURSE_SOURCE.exists() and TIMECOURSE_EXPORT.exists():
        if TIMECOURSE_SOURCE.stat().st_mtime > TIMECOURSE_EXPORT.stat().st_mtime:
            issues.append(
                "Timecourse export is older than its source H5AD; regenerate the export"
            )

    print("\nOptional analysis coverage")
    for label, path in OPTIONAL_ASSETS.items():
        if path.is_dir():
            count = len(list(path.glob("*.csv"))) + len(list(path.glob("*.png")))
            print(f"{'✓' if count else '–'} {label}: {count:,} files")
        elif path.exists():
            print(f"✓ {label}")
        else:
            print(f"– {label}: not exported")

    if issues:
        print("\nRelease-blocking issues")
        for issue in issues:
            print(f"✗ {issue}")
        return 1

    print("\n✓ EBs portal data is ready")
    return 0


if __name__ == "__main__":
    sys.exit(main())
