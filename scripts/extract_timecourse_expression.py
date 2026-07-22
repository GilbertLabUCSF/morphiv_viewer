#!/usr/bin/env python3
"""
Extract per-gene expression statistics across timepoints from the full merged
timecourse matrix.

This script reads the full merged timecourse h5ad file and extracts:
- Mean log1p-normalized counts per 10,000 (CP10K) per gene per day
- Percentage of cells with observed expression (>0)
- Number of cells per day

Output: data_extracted/timecourse_expression.parquet
"""

from pathlib import Path
from datetime import datetime, timezone
import json
import warnings

import numpy as np
import pandas as pd

# Suppress warnings during import
warnings.filterwarnings("ignore", category=FutureWarning)

try:
    import scanpy as sc
except ImportError:
    print("ERROR: scanpy is required for this script.")
    print("Install with: pip install scanpy")
    exit(1)

# Paths
SCRIPT_DIR = Path(__file__).parent
PROJECT_ROOT = SCRIPT_DIR.parent
DATA_SYMLINK = PROJECT_ROOT / "data"
DATA_EXTRACTED = PROJECT_ROOT / "data_extracted"

# Prefer the full-gene merged object. The scVI object contains only selected
# HVGs/markers and would omit most screened transcription factors.
H5AD_PATHS = [
    DATA_SYMLINK / "single_cell_timecourse/results/processed/timecourse_merged.h5ad",
    DATA_SYMLINK / "single_cell_timecourse/results_latest/processed/timecourse_merged.h5ad",
]

OUTPUT_PATH = DATA_EXTRACTED / "timecourse_expression.parquet"
MANIFEST_PATH = DATA_EXTRACTED / "timecourse_manifest.json"

# Column that contains day/timepoint information
DAY_COLUMN = "day"  # Adjust if different in the actual data


def find_h5ad():
    """Find the timecourse h5ad file."""
    for path in H5AD_PATHS:
        if path.exists():
            return path
    return None


def extract_expression_stats(adata, day_column: str) -> pd.DataFrame:
    """
    Extract per-gene expression statistics grouped by day.

    Returns DataFrame with columns:
    - gene: gene symbol
    - day: timepoint
    - mean_expression: mean log1p-normalized CP10K expression
    - pct_detected: percentage of cells with observed expression > 0
    - n_cells: number of cells in this day
    """
    print(f"  Extracting expression stats grouped by '{day_column}'...")

    # Check if day column exists
    if day_column not in adata.obs.columns:
        # Try to find alternative column names
        possible_cols = [c for c in adata.obs.columns if "day" in c.lower() or "time" in c.lower()]
        if possible_cols:
            day_column = possible_cols[0]
            print(f"  Using alternative column: {day_column}")
        else:
            print(f"  Available columns: {list(adata.obs.columns[:20])}...")
            raise ValueError(f"Cannot find day/timepoint column. Tried: {DAY_COLUMN}")

    # Get unique days
    days = adata.obs[day_column].unique()
    print(f"  Found {len(days)} timepoints: {sorted(days)}")

    # The merged pipeline stores log1p(CP10K) in X. Zeros remain zeros, so the
    # same sparse matrix supports both a mean and an observed-detection rate.
    X = adata.X
    print("  Using X: log1p-normalized counts per 10,000")

    # Get gene names
    genes = adata.var_names.tolist()

    results = []

    for day in sorted(days):
        mask = (adata.obs[day_column] == day).to_numpy()
        indices = np.flatnonzero(mask)
        n_cells = len(indices)

        if n_cells == 0:
            continue

        # Subset data for this day
        X_day = X[indices, :]

        # Keep the operation sparse; densifying a full day can require several GB.
        mean_expr = np.asarray(X_day.mean(axis=0)).ravel()
        pct_detected = np.asarray((X_day > 0).mean(axis=0)).ravel() * 100

        # Create rows for this day
        for i, gene in enumerate(genes):
            results.append({
                "gene": gene,
                "day": day,
                "mean_expression": float(mean_expr[i]),
                "pct_detected": float(pct_detected[i]),
                # Backward-compatible alias for downstream downloads.
                "pct_expressing": float(pct_detected[i]),
                "n_cells": int(n_cells),
                "expression_scale": "log1p_cp10k",
            })

        print(f"    Day {day}: {n_cells} cells processed")

    df = pd.DataFrame(results)
    return df


def main():
    print("=" * 60)
    print("Extracting Timecourse Expression Data")
    print("=" * 60)
    print()

    # Check data symlink
    if not DATA_SYMLINK.exists():
        print("ERROR: Data symlink does not exist.")
        print(f"  Expected: {DATA_SYMLINK}")
        print("  Run: ln -s ../morphic_pub_refactor data")
        exit(1)

    # Find h5ad file
    h5ad_path = find_h5ad()
    if h5ad_path is None:
        print("ERROR: Cannot find timecourse h5ad file.")
        print("  Searched locations:")
        for path in H5AD_PATHS:
            print(f"    - {path}")
        exit(1)

    print(f"[1/4] Found h5ad: {h5ad_path}")
    print(f"  File size: {h5ad_path.stat().st_size / 1e9:.2f} GB")
    print()

    # Load h5ad
    print("[2/4] Loading h5ad file (this may take a few minutes)...")
    adata = sc.read_h5ad(h5ad_path, backed="r")
    print(f"  Loaded: {adata.n_obs} cells x {adata.n_vars} genes")
    print(f"  Obs columns: {list(adata.obs.columns[:10])}...")
    print()

    # Extract expression stats
    print("[3/4] Calculating expression statistics...")
    df = extract_expression_stats(adata, DAY_COLUMN)
    print(f"  Generated {len(df)} rows")
    print()

    # Save output
    print("[4/4] Saving to parquet...")
    DATA_EXTRACTED.mkdir(exist_ok=True)
    temp_path = OUTPUT_PATH.with_suffix(".parquet.tmp")
    df.to_parquet(temp_path, index=False)
    temp_path.replace(OUTPUT_PATH)
    manifest = {
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "source": str(h5ad_path.relative_to(PROJECT_ROOT)),
        "source_modified_at": datetime.fromtimestamp(
            h5ad_path.stat().st_mtime, timezone.utc
        ).isoformat(),
        "n_cells": int(adata.n_obs),
        "n_genes": int(adata.n_vars),
        "timepoints": [str(value) for value in sorted(adata.obs[DAY_COLUMN].unique())],
        "expression_scale": "log1p-normalized counts per 10,000",
        "detection_definition": "fraction of cells with observed expression > 0",
    }
    MANIFEST_PATH.write_text(json.dumps(manifest, indent=2) + "\n")
    adata.file.close()
    print(f"  Saved: {OUTPUT_PATH}")
    print(f"  Size: {OUTPUT_PATH.stat().st_size / 1e6:.2f} MB")
    print()

    print("=" * 60)
    print("✓ Extraction complete!")
    print("=" * 60)

    # Show sample output
    print()
    print("Sample output:")
    print(df.head(10).to_string())


if __name__ == "__main__":
    main()
