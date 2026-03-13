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
    DEG_TABLES_EBS,
    DEG_TABLES_IPSC,
    DOSE_RESPONSE_EBS,
    DOSE_RESPONSE_IPSC,
    LINEAGE_DE_EBS,
    LINEAGE_DE_IPSC,
    PATHWAY_ENRICHMENT_EBS,
    PATHWAY_ENRICHMENT_IPSC,
    RESULTS_EBS,
    RESULTS_IPSC,
    TF_PERTURBSEQ,
    TF_SIMILARITY_EBS,
    TF_SIMILARITY_IPSC,
    TRANSCRIPTOME_EDIST_EBS,
    TRANSCRIPTOME_EDIST_IPSC,
    VIABILITY_EBS,
    VIABILITY_IPSC,
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
        raise FileNotFoundError("No compositional analysis files found")

    df = pd.concat(dfs, ignore_index=True)

    # Extract gene symbol
    df["gene"] = df["perturbation"].str.split("_").str[0]

    return df


def load_deg_table(perturbation: str, condition: str = None) -> Optional[pd.DataFrame]:
    """
    Load DEG table for a specific perturbation.

    Args:
        perturbation: Perturbation ID (e.g., "SOX2_P1P2")
        condition: Optional condition to search ("EBs" or "iPSC").
                   If None, searches both and returns first found.

    Returns DataFrame with columns:
    - gene, baseMean, log2FoldChange, lfcSE, pvalue, padj
    - significant, is_deg

    Returns None if no DEG data exists for this perturbation.
    """
    search_dirs = []
    if condition == "EBs":
        search_dirs = [DEG_TABLES_EBS]
    elif condition == "iPSC":
        search_dirs = [DEG_TABLES_IPSC]
    else:
        search_dirs = [DEG_TABLES_EBS, DEG_TABLES_IPSC]

    for deg_dir in search_dirs:
        if not deg_dir.exists():
            continue

        deg_path = deg_dir / f"{perturbation}_deg.csv"
        if deg_path.exists():
            return pd.read_csv(deg_path)

        # Try glob match
        possible_files = list(deg_dir.glob(f"{perturbation}*_deg.csv"))
        if possible_files:
            return pd.read_csv(possible_files[0])

    return None


@st.cache_data(ttl=3600)
def load_dose_response(condition: str) -> Optional[pd.DataFrame]:
    """
    Load dose response data for a condition.

    Returns DataFrame with columns:
    - perturbation, knockdown_pct, lineage_effect, e_distance
    """
    dr_dir = DOSE_RESPONSE_EBS if condition == "EBs" else DOSE_RESPONSE_IPSC
    path = dr_dir / f"{condition}_dose_response.csv"

    if not path.exists():
        return None

    return pd.read_csv(path)


@st.cache_data(ttl=3600)
def load_pathway_enrichment(condition: str) -> Optional[pd.DataFrame]:
    """
    Load pathway enrichment results for a condition.

    Returns DataFrame with columns:
    - perturbation, term, gene_set_library, p_value, adjusted_p_value,
      overlap, genes, Odds Ratio, Combined Score
    """
    pe_dir = PATHWAY_ENRICHMENT_EBS if condition == "EBs" else PATHWAY_ENRICHMENT_IPSC
    path = pe_dir / f"{condition}_enrichment_results.csv"

    if not path.exists():
        return None

    return pd.read_csv(path)


@st.cache_data(ttl=3600)
def load_transcriptome_edist(condition: str) -> Optional[pd.DataFrame]:
    """
    Load transcriptome E-distance data for a condition.

    Returns DataFrame with columns:
    - perturbation, n_cells, edist_observed, edist_null_mean, edist_null_std,
      p_value, q_value
    """
    ed_dir = TRANSCRIPTOME_EDIST_EBS if condition == "EBs" else TRANSCRIPTOME_EDIST_IPSC
    path = ed_dir / f"{condition}_tedist_merged.csv"

    if not path.exists():
        return None

    return pd.read_csv(path)


@st.cache_data(ttl=3600)
def load_tf_clusters(condition: str) -> Optional[pd.DataFrame]:
    """
    Load TF cluster assignments for a condition.

    Returns DataFrame with columns:
    - perturbation, cluster
    """
    tf_dir = TF_SIMILARITY_EBS if condition == "EBs" else TF_SIMILARITY_IPSC
    path = tf_dir / f"{condition}_tf_clusters.csv"

    if not path.exists():
        return None

    return pd.read_csv(path)


def load_lineage_de(perturbation: str, lineage: str, condition: str) -> Optional[pd.DataFrame]:
    """
    Load lineage-specific DE results for a perturbation.

    Args:
        perturbation: Perturbation ID (e.g., "SOX2_P1P2")
        lineage: Lineage name (e.g., "Neural_Ectoderm")
        condition: "EBs" or "iPSC"

    Returns DataFrame with columns:
    - gene (ENSG IDs), baseMean, log2fc, lfcSE, stat, pval, padj, is_deg
    """
    de_dir = LINEAGE_DE_EBS if condition == "EBs" else LINEAGE_DE_IPSC
    path = de_dir / f"{perturbation}_{lineage}_DE.csv"

    if not path.exists():
        return None

    return pd.read_csv(path)


def load_double_diff(perturbation: str, lineage1: str, lineage2: str, condition: str) -> Optional[pd.DataFrame]:
    """
    Load double differential results between two lineages.

    Returns DataFrame with columns:
    - double_diff, abs_double_diff, is_significant
    """
    de_dir = LINEAGE_DE_EBS if condition == "EBs" else LINEAGE_DE_IPSC
    path = de_dir / f"{perturbation}_{lineage1}_vs_{lineage2}_double_diff.csv"

    if not path.exists():
        return None

    return pd.read_csv(path)


def get_available_lineage_de(perturbation: str, condition: str) -> list[str]:
    """
    Get list of lineages with DE data for a perturbation.

    Checks which {perturbation}_{lineage}_DE.csv files exist.
    """
    de_dir = LINEAGE_DE_EBS if condition == "EBs" else LINEAGE_DE_IPSC

    if not de_dir.exists():
        return []

    lineages = []
    for f in de_dir.glob(f"{perturbation}_*_DE.csv"):
        # Extract lineage name: {perturbation}_{lineage}_DE.csv
        name = f.stem  # e.g., "SOX2_P1P2_Neural_Ectoderm_DE"
        # Remove perturbation prefix and _DE suffix
        suffix = name[len(perturbation) + 1:]  # e.g., "Neural_Ectoderm_DE"
        if suffix.endswith("_DE"):
            lineage = suffix[:-3]  # e.g., "Neural_Ectoderm"
            lineages.append(lineage)

    return sorted(lineages)


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
            "Timecourse expression not extracted. "
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
