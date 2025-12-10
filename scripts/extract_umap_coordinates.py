#!/usr/bin/env python3
"""
Extract UMAP coordinates and cell metadata for interactive plotting.

This extracts a lightweight version of the UMAP data that can be loaded
quickly for dynamic highlighting of perturbations.

Output: data_extracted/umap_coordinates_{condition}.parquet
"""

from pathlib import Path
import warnings
import sys

import numpy as np
import pandas as pd

warnings.filterwarnings("ignore", category=FutureWarning)

try:
    import scanpy as sc
except ImportError:
    print("ERROR: scanpy is required. Install with: pip install scanpy")
    sys.exit(1)

# Paths
SCRIPT_DIR = Path(__file__).parent
PROJECT_ROOT = SCRIPT_DIR.parent
DATA_SYMLINK = PROJECT_ROOT / "data"
DATA_EXTRACTED = PROJECT_ROOT / "data_extracted"

# Use the compact version which should be smaller
H5AD_PATHS = {
    "EBs": DATA_SYMLINK / "tf_perturbseq/results/EBs/cell_type_scored_scvi_compact.h5ad",
    "iPSC": DATA_SYMLINK / "tf_perturbseq/results/iPSC/cell_type_scored_scvi_compact.h5ad",
}

# Fallback to regular if compact doesn't exist
H5AD_FALLBACK = {
    "EBs": DATA_SYMLINK / "tf_perturbseq/results/EBs/cell_type_scored_scvi.h5ad",
    "iPSC": DATA_SYMLINK / "tf_perturbseq/results/iPSC/cell_type_scored_scvi.h5ad",
}


def extract_umap_for_condition(condition: str):
    """Extract UMAP coordinates and metadata for a condition."""

    h5ad_path = H5AD_PATHS.get(condition)
    if not h5ad_path or not h5ad_path.exists():
        h5ad_path = H5AD_FALLBACK.get(condition)

    if not h5ad_path or not h5ad_path.exists():
        print(f"  ✗ No h5ad found for {condition}")
        return None

    print(f"  Loading {h5ad_path.name}...")
    print(f"    Size: {h5ad_path.stat().st_size / 1e9:.2f} GB")

    # Load with backed mode to reduce memory
    adata = sc.read_h5ad(h5ad_path, backed='r')

    print(f"    Cells: {adata.n_obs:,}")

    # Check for UMAP coordinates
    umap_key = None
    for key in ['X_umap', 'X_umap_scvi', 'X_umap_scanvi']:
        if key in adata.obsm:
            umap_key = key
            break

    if umap_key is None:
        print(f"  ✗ No UMAP coordinates found in obsm: {list(adata.obsm.keys())}")
        return None

    print(f"    Using UMAP: {umap_key}")

    # Extract UMAP coordinates
    umap_coords = adata.obsm[umap_key]

    # Build dataframe with essential columns
    df = pd.DataFrame({
        'umap_1': umap_coords[:, 0],
        'umap_2': umap_coords[:, 1],
    })

    # Add perturbation info
    for col in ['perturbation', 'gene', 'target_gene', 'guide']:
        if col in adata.obs.columns:
            df[col] = adata.obs[col].values

    # Add cell type info
    for col in ['cell_type', 'predicted_cell_type', 'scanvi_pred', 'dominant_cell_type']:
        if col in adata.obs.columns:
            df['cell_type'] = adata.obs[col].values
            break

    # Add lineage scores if available
    lineage_cols = [c for c in adata.obs.columns if c.endswith('_score') or c in
                    ['Amnion', 'Epiblast', 'Formative_Epiblast', 'Neural_Ectoderm',
                     'Non_neural_Ectoderm', 'Trophoblast_Like']]
    for col in lineage_cols[:6]:  # Limit to 6 lineages
        if col in adata.obs.columns:
            df[col] = adata.obs[col].values

    # Sample if too many cells (for performance)
    max_cells = 100000
    if len(df) > max_cells:
        print(f"    Sampling {max_cells:,} cells from {len(df):,}")
        df = df.sample(n=max_cells, random_state=42)

    return df


def main():
    print("=" * 60)
    print("Extracting UMAP Coordinates")
    print("=" * 60)
    print()

    if not DATA_SYMLINK.exists():
        print("ERROR: Data symlink does not exist")
        sys.exit(1)

    DATA_EXTRACTED.mkdir(exist_ok=True)

    for condition in ["EBs", "iPSC"]:
        print(f"[{condition}]")

        df = extract_umap_for_condition(condition)

        if df is not None:
            output_path = DATA_EXTRACTED / f"umap_coordinates_{condition}.parquet"
            df.to_parquet(output_path, index=False)
            print(f"    Saved: {output_path}")
            print(f"    Size: {output_path.stat().st_size / 1e6:.2f} MB")
            print(f"    Rows: {len(df):,}")
        print()

    print("=" * 60)
    print("Done!")
    print("=" * 60)


if __name__ == "__main__":
    main()
