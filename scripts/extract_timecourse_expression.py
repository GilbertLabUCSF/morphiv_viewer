#!/usr/bin/env python3
"""
Extract per-gene expression statistics across timepoints from timecourse h5ad.

This script reads the timecourse scVI-processed h5ad file and extracts:
- Mean expression per gene per day
- Percentage of cells expressing at multiple thresholds:
  - pct_expressing: >0 (any expression)
  - pct_expr_gt1: >1 TPM (low threshold)
  - pct_expr_gt5: >5 TPM (moderate threshold)
- Number of cells per day

Output: data_extracted/timecourse_expression.parquet
"""

from pathlib import Path
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

# Input h5ad - try multiple possible locations
H5AD_PATHS = [
    DATA_SYMLINK / "single_cell_timecourse/results_latest/processed/timecourse_scvi.h5ad",
    DATA_SYMLINK / "single_cell_timecourse/results/processed/timecourse_scvi.h5ad",
]

OUTPUT_PATH = DATA_EXTRACTED / "timecourse_expression.parquet"

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
    - mean_expression: mean expression (from normalized layer or X)
    - pct_expressing: percentage of cells with expression > 0
    - pct_expr_gt1: percentage of cells with expression > 1 TPM
    - pct_expr_gt5: percentage of cells with expression > 5 TPM
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

    # Use normalized counts if available, otherwise X
    if "X_normalized" in adata.layers:
        X = adata.layers["X_normalized"]
        print("  Using X_normalized layer")
    elif "counts" in adata.layers:
        # If we have counts, use X (should be normalized)
        X = adata.X
        print("  Using X matrix")
    else:
        X = adata.X
        print("  Using X matrix (default)")

    # Get gene names
    genes = adata.var_names.tolist()

    results = []

    for day in sorted(days):
        mask = adata.obs[day_column] == day
        n_cells = mask.sum()

        if n_cells == 0:
            continue

        # Subset data for this day
        X_day = X[mask, :]

        # Handle sparse matrix
        if hasattr(X_day, "toarray"):
            X_day_dense = X_day.toarray()
        else:
            X_day_dense = np.array(X_day)

        # Calculate statistics
        mean_expr = np.mean(X_day_dense, axis=0)
        pct_expr = np.mean(X_day_dense > 0, axis=0) * 100
        pct_expr_gt1 = np.mean(X_day_dense > 1, axis=0) * 100
        pct_expr_gt5 = np.mean(X_day_dense > 5, axis=0) * 100

        # Create rows for this day
        for i, gene in enumerate(genes):
            results.append({
                "gene": gene,
                "day": day,
                "mean_expression": float(mean_expr[i]),
                "pct_expressing": float(pct_expr[i]),
                "pct_expr_gt1": float(pct_expr_gt1[i]),
                "pct_expr_gt5": float(pct_expr_gt5[i]),
                "n_cells": int(n_cells),
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
    adata = sc.read_h5ad(h5ad_path)
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
    df.to_parquet(OUTPUT_PATH, index=False)
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
