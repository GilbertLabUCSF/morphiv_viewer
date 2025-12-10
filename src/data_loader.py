"""
Data loading functions with Streamlit caching.

All loaders check for data existence and raise clear errors if missing.
"""

from pathlib import Path
from typing import Optional

import pandas as pd
import streamlit as st

from .config import (
    DATA_EXTRACTED,
    DEG_TABLES,
    RESULTS_EBS,
    RESULTS_IPSC,
    VIABILITY_EBS,
    VIABILITY_IPSC,
    TF_PERTURBSEQ,
)


@st.cache_data(ttl=3600)
def load_lineage_analysis() -> pd.DataFrame:
    """
    Load lineage analysis results for both conditions.

    Returns DataFrame with columns including:
    - perturbation, condition, lineage, n_cells
    - observed_*, z_score_*, p_value_* for various metrics
    - effect_size_interpretation, max_abs_z_score
    """
    ebs_path = RESULTS_EBS / "lineage_analysis" / "EBs_lineage_analysis_merged.csv"
    ipsc_path = RESULTS_IPSC / "lineage_analysis" / "iPSC_lineage_analysis_merged.csv"

    dfs = []
    missing = []

    if ebs_path.exists():
        ebs = pd.read_csv(ebs_path)
        ebs["condition"] = "EBs"
        dfs.append(ebs)
    else:
        missing.append(f"EBs ({ebs_path})")

    if ipsc_path.exists():
        ipsc = pd.read_csv(ipsc_path)
        ipsc["condition"] = "iPSC"
        dfs.append(ipsc)
    else:
        missing.append(f"iPSC ({ipsc_path})")

    if not dfs:
        raise FileNotFoundError(
            "Lineage analysis files not found. "
            f"Expected: {ebs_path} or {ipsc_path}"
        )

    df = pd.concat(dfs, ignore_index=True)

    if missing:
        st.warning(f"Lineage analysis missing for: {', '.join(missing)}. Showing available data only.")

    # Extract gene symbol from perturbation (e.g., "SOX2_P1P2" -> "SOX2")
    df["gene"] = df["perturbation"].str.split("_").str[0]

    return df


@st.cache_data(ttl=3600)
def load_knockdown_efficiency() -> pd.DataFrame:
    """
    Load knockdown efficiency data for both conditions.

    Returns DataFrame with columns:
    - perturbation, target_gene, condition
    - log2_fold_change, knockdown_pct, knockdown_category
    - dominant_cell_type, cell_type_purity, n_cells
    """
    ebs_path = RESULTS_EBS / "knockdown_efficiency" / "knockdown_efficiency_all_genes.csv"
    ipsc_path = RESULTS_IPSC / "knockdown_efficiency" / "knockdown_efficiency_all_genes.csv"

    dfs = []
    missing = []

    if ebs_path.exists():
        ebs = pd.read_csv(ebs_path)
        ebs["condition"] = "EBs"
        dfs.append(ebs)
    else:
        missing.append(f"EBs ({ebs_path})")

    if ipsc_path.exists():
        ipsc = pd.read_csv(ipsc_path)
        ipsc["condition"] = "iPSC"
        dfs.append(ipsc)
    else:
        missing.append(f"iPSC ({ipsc_path})")

    if not dfs:
        raise FileNotFoundError(
            "Knockdown efficiency data not found. "
            f"Expected: {ebs_path} or {ipsc_path}"
        )

    df = pd.concat(dfs, ignore_index=True)

    if missing:
        st.warning(f"Knockdown efficiency missing for: {', '.join(missing)}. Showing available data only.")
    return df


@st.cache_data(ttl=3600)
def load_viability() -> pd.DataFrame:
    """
    Load viability/fitness scores for both conditions.

    Returns DataFrame with columns:
    - gene, condition
    - median_lfc, mean_lfc, stouffer_z, stouffer_p
    - n_depleted, frac_depleted
    """
    ebs_path = VIABILITY_EBS / "viability_scores_gene_level.csv"
    ipsc_path = VIABILITY_IPSC / "viability_scores_gene_level.csv"

    dfs = []
    missing = []

    if ebs_path.exists():
        ebs = pd.read_csv(ebs_path)
        ebs["condition"] = "EBs"
        dfs.append(ebs)
    else:
        missing.append(f"EBs ({ebs_path})")

    if ipsc_path.exists():
        ipsc = pd.read_csv(ipsc_path)
        ipsc["condition"] = "iPSC"
        dfs.append(ipsc)
    else:
        missing.append(f"iPSC ({ipsc_path})")

    if not dfs:
        raise FileNotFoundError(
            "Viability data not found. "
            f"Expected: {ebs_path} or {ipsc_path}"
        )

    df = pd.concat(dfs, ignore_index=True)

    if missing:
        st.warning(f"Viability data missing for: {', '.join(missing)}. Showing available data only.")
    return df


@st.cache_data(ttl=3600)
def load_compositional() -> pd.DataFrame:
    """
    Load compositional analysis (chi-square) results for both conditions.

    Returns DataFrame with columns:
    - perturbation, condition
    - chi2_statistic, p_value, q_value
    - residual_* for each lineage
    """
    ebs_path = RESULTS_EBS / "compositional" / "EBs_significant_hits.csv"
    ipsc_path = RESULTS_IPSC / "compositional" / "iPSC_significant_hits.csv"

    dfs = []

    if ebs_path.exists():
        ebs = pd.read_csv(ebs_path)
        ebs["condition"] = "EBs"
        dfs.append(ebs)

    if ipsc_path.exists():
        ipsc = pd.read_csv(ipsc_path)
        ipsc["condition"] = "iPSC"
        dfs.append(ipsc)

    if not dfs:
        raise FileNotFoundError("NEEDS IMPLEMENTATION: No compositional analysis files found")

    df = pd.concat(dfs, ignore_index=True)

    # Extract gene symbol
    df["gene"] = df["perturbation"].str.split("_").str[0]

    return df


def load_deg_table(perturbation: str) -> Optional[pd.DataFrame]:
    """
    Load DEG table for a specific perturbation.

    Returns DataFrame with columns:
    - gene, baseMean, log2FoldChange, lfcSE, pvalue, padj
    - significant, is_deg

    Returns None if no DEG data exists for this perturbation.
    """
    # Try exact match first
    deg_path = DEG_TABLES / f"{perturbation}_deg.csv"

    if not deg_path.exists():
        # Try without the isoform suffix (e.g., "SOX2_P1P2" -> check for "SOX2_P1P2_deg.csv")
        possible_files = list(DEG_TABLES.glob(f"{perturbation}*_deg.csv"))
        if possible_files:
            deg_path = possible_files[0]
        else:
            return None

    if not deg_path.exists():
        return None

    df = pd.read_csv(deg_path)
    return df


@st.cache_data(ttl=3600)
def load_timecourse_expression() -> pd.DataFrame:
    """
    Load pre-extracted timecourse expression data.

    Returns DataFrame with columns:
    - gene, day, mean_expression, pct_expressing, n_cells

    Raises FileNotFoundError if extraction script hasn't been run.
    """
    path = DATA_EXTRACTED / "timecourse_expression.parquet"

    if not path.exists():
        raise FileNotFoundError(
            "NEEDS IMPLEMENTATION: Timecourse expression not extracted. "
            "Run: python scripts/extract_timecourse_expression.py"
        )

    return pd.read_parquet(path)


@st.cache_data(ttl=3600)
def get_gene_list() -> list[str]:
    """
    Get sorted list of all unique gene symbols from the screen.
    """
    # Use resolved targets from EBs (should be same genes)
    ebs_targets = RESULTS_EBS / "resolved_targets.csv"
    ipsc_targets = RESULTS_IPSC / "resolved_targets.csv"

    genes = set()

    if ebs_targets.exists():
        df = pd.read_csv(ebs_targets)
        # The file might have different column names - try common ones
        for col in ["target_gene", "gene", "gene_symbol", "target"]:
            if col in df.columns:
                genes.update(df[col].dropna().unique())
                break

    if ipsc_targets.exists():
        df = pd.read_csv(ipsc_targets)
        for col in ["target_gene", "gene", "gene_symbol", "target"]:
            if col in df.columns:
                genes.update(df[col].dropna().unique())
                break

    # If no targets file, extract from lineage analysis
    if not genes:
        try:
            la = load_lineage_analysis()
            genes = set(la["gene"].unique())
        except FileNotFoundError:
            pass

    return sorted(genes)


@st.cache_data(ttl=3600)
def get_perturbations_for_gene(gene: str) -> dict:
    """
    Get all perturbation IDs for a given gene.

    Returns dict with:
    - perturbations: list of perturbation IDs
    - conditions: dict mapping perturbation to list of conditions tested
    """
    la = load_lineage_analysis()
    gene_data = la[la["gene"] == gene]

    perturbations = gene_data["perturbation"].unique().tolist()

    conditions = {}
    for pert in perturbations:
        pert_data = gene_data[gene_data["perturbation"] == pert]
        conditions[pert] = pert_data["condition"].unique().tolist()

    return {
        "perturbations": perturbations,
        "conditions": conditions,
    }


@st.cache_data(ttl=3600)
def get_summary_stats() -> dict:
    """
    Get summary statistics for the portal landing page.
    Cell counts are from the h5ad files (total cells in dataset).
    """
    stats = {}

    try:
        la = load_lineage_analysis()
        stats["n_genes"] = la["gene"].nunique()
        stats["n_perturbations"] = la["perturbation"].nunique()

        # Cell counts from h5ad files (accurate totals)
        stats["n_cells_ebs"] = 952_347
        stats["n_cells_ipsc"] = 574_589
        stats["n_cells_timecourse"] = 91_185
        stats["n_cells_total"] = stats["n_cells_ipsc"] + stats["n_cells_ebs"]

    except FileNotFoundError as e:
        stats["error"] = str(e)

    return stats


def load_marker_counts(perturbation: str, condition: str) -> Optional[pd.DataFrame]:
    """
    Load marker expression counts for a specific perturbation.

    Used for generating dynamic dotplots.

    Returns DataFrame with columns:
    - gene, raw_counts, cpm, n_cells, condition

    Returns None if no data exists for this perturbation.
    """
    marker_dir = TF_PERTURBSEQ / "results" / condition / "marker_pseudobulk" / "per_guide"
    marker_path = marker_dir / f"{perturbation}_marker_counts.csv"

    if not marker_path.exists():
        return None

    df = pd.read_csv(marker_path)
    return df


@st.cache_data(ttl=86400, show_spinner="Loading UMAP coordinates...")
def load_umap_coordinates(condition: str, max_cells: int = 50000) -> Optional[pd.DataFrame]:
    """
    Load UMAP coordinates for a condition.

    First tries pre-extracted parquet (fast), then falls back to h5ad (slow but works).
    Uses subsampling to keep memory/speed reasonable.

    Returns DataFrame with columns:
    - umap_1, umap_2, perturbation, cell_type, etc.

    Returns None if no data available.
    """
    # Try pre-extracted first (fast path)
    path = DATA_EXTRACTED / f"umap_coordinates_{condition}.parquet"
    if path.exists():
        return pd.read_parquet(path)

    # Fall back to loading from h5ad (slow but works)
    h5ad_paths = {
        "EBs": TF_PERTURBSEQ / "results" / "EBs" / "cell_type_scored_scvi_compact.h5ad",
        "iPSC": TF_PERTURBSEQ / "results" / "iPSC" / "cell_type_scored_scvi_compact.h5ad",
    }

    h5ad_path = h5ad_paths.get(condition)
    if not h5ad_path or not h5ad_path.exists():
        # Try non-compact version
        h5ad_path = TF_PERTURBSEQ / "results" / condition / "cell_type_scored_scvi.h5ad"

    if not h5ad_path or not h5ad_path.exists():
        return None

    try:
        import scanpy as sc
        import numpy as np

        # Load with backed mode to reduce memory
        adata = sc.read_h5ad(h5ad_path, backed='r')

        # Find UMAP key
        umap_key = None
        for key in ['X_umap', 'X_umap_scvi', 'X_umap_scanvi']:
            if key in adata.obsm:
                umap_key = key
                break

        if umap_key is None:
            return None

        # Subsample if too many cells
        n_cells = adata.n_obs
        if n_cells > max_cells:
            indices = np.random.choice(n_cells, max_cells, replace=False)
            indices = np.sort(indices)
        else:
            indices = np.arange(n_cells)

        # Extract data
        umap_coords = adata.obsm[umap_key][indices]

        df = pd.DataFrame({
            'umap_1': umap_coords[:, 0],
            'umap_2': umap_coords[:, 1],
        })

        # Add perturbation info
        for col in ['perturbation', 'gene', 'target_gene', 'guide']:
            if col in adata.obs.columns:
                df[col] = adata.obs[col].values[indices]

        # Add cell type - try direct column first
        cell_type_found = False
        for col in ['cell_type', 'predicted_cell_type', 'scanvi_predictions', 'scanvi_pred', 'dominant_cell_type']:
            if col in adata.obs.columns:
                df['cell_type'] = adata.obs[col].values[indices]
                cell_type_found = True
                break

        # If no cell type, derive from lineage scores
        if not cell_type_found:
            score_cols = [c for c in adata.obs.columns if c.endswith('_score') and not c.endswith('_score_raw')]
            if score_cols:
                import numpy as np
                scores = np.column_stack([adata.obs[c].values[indices] for c in score_cols])
                lineage_names = [c.replace('_score', '') for c in score_cols]
                dominant_idx = np.argmax(scores, axis=1)
                df['cell_type'] = [lineage_names[i] for i in dominant_idx]

        return df

    except ImportError:
        return None
    except Exception:
        return None


@st.cache_data(ttl=3600)
def load_all_marker_counts(condition: str) -> Optional[pd.DataFrame]:
    """
    Load all marker counts for a condition into a combined DataFrame.

    Used for efficient dotplot generation across perturbations.
    """
    marker_dir = TF_PERTURBSEQ / "results" / condition / "marker_pseudobulk" / "per_guide"

    if not marker_dir.exists():
        return None

    all_dfs = []
    for f in marker_dir.glob("*_marker_counts.csv"):
        df = pd.read_csv(f)
        # Extract perturbation from filename
        pert = f.stem.replace("_marker_counts", "")
        df["perturbation"] = pert
        all_dfs.append(df)

    if not all_dfs:
        return None

    return pd.concat(all_dfs, ignore_index=True)
