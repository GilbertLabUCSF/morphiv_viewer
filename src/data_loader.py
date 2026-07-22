"""Cached data access for the EBs release browser."""

from __future__ import annotations

from datetime import datetime
import json
from typing import Optional

import pandas as pd
import streamlit as st

from .config import (
    DATA_EXTRACTED,
    DEG_TABLES_EBS,
    DOSE_RESPONSE_EBS,
    LINEAGE_DE_EBS,
    PATHWAY_ENRICHMENT_EBS,
    RESULTS_EBS,
    TF_SIMILARITY_EBS,
    TRANSCRIPTOME_EDIST_EBS,
    VIABILITY_EBS,
)


def _attach_target_genes(df: pd.DataFrame) -> pd.DataFrame:
    """Attach the screen's canonical target symbol to perturbation-level rows."""
    if "perturbation" not in df.columns:
        return df
    gene_map = (
        load_lineage_analysis()
        .drop_duplicates("perturbation")
        .set_index("perturbation")["gene"]
    )
    result = df.copy()
    result["gene"] = result["perturbation"].map(gene_map)
    result["gene"] = result["gene"].fillna(
        result["perturbation"].str.split("_").str[0]
    )
    return result


@st.cache_data(ttl=3600)
def load_lineage_analysis() -> pd.DataFrame:
    """Load the primary EBs lineage-effect table with canonical target symbols."""
    path = RESULTS_EBS / "lineage_analysis" / "EBs_lineage_analysis_merged.csv"
    if not path.exists():
        raise FileNotFoundError(f"EB lineage analysis not found: {path}")

    result = pd.read_csv(path)
    result["condition"] = "EBs"
    result["source_gene"] = result["perturbation"].str.split("_").str[0]

    target_path = (
        RESULTS_EBS
        / "knockdown_efficiency"
        / "knockdown_efficiency_all_genes.csv"
    )
    if target_path.exists():
        targets = pd.read_csv(
            target_path, usecols=["perturbation", "target_gene"]
        ).drop_duplicates("perturbation")
        result = result.merge(targets, on="perturbation", how="left")
        result["gene"] = result["target_gene"].fillna(result["source_gene"])
    else:
        result["gene"] = result["source_gene"]
    return result


@st.cache_data(ttl=3600)
def load_knockdown_efficiency() -> pd.DataFrame:
    """Load EBs knockdown-efficiency estimates."""
    path = (
        RESULTS_EBS
        / "knockdown_efficiency"
        / "knockdown_efficiency_all_genes.csv"
    )
    if not path.exists():
        raise FileNotFoundError(f"EB knockdown efficiency not found: {path}")
    result = pd.read_csv(path)
    result["condition"] = "EBs"
    return result


@st.cache_data(ttl=3600)
def load_viability() -> pd.DataFrame:
    """Load EBs gene-level viability/fitness scores."""
    path = VIABILITY_EBS / "viability_scores_gene_level.csv"
    if not path.exists():
        raise FileNotFoundError(f"EB viability data not found: {path}")
    result = pd.read_csv(path)
    result["condition"] = "EBs"
    return result


@st.cache_data(ttl=3600)
def load_compositional() -> pd.DataFrame:
    """Load significant EBs compositional probability-shift hits."""
    path = RESULTS_EBS / "compositional" / "EBs_significant_hits.csv"
    if not path.exists():
        raise FileNotFoundError(f"EB compositional results not found: {path}")
    result = _attach_target_genes(pd.read_csv(path))
    result["condition"] = "EBs"
    return result


@st.cache_data(ttl=3600)
def load_probability_shifts() -> pd.DataFrame:
    """Load EBs k-nearest-neighbor lineage probability shifts."""
    path = RESULTS_EBS / "knn_probability_shifts" / "EBs_knn_prob_shifts.csv"
    if not path.exists():
        raise FileNotFoundError(f"EB probability shifts not found: {path}")
    result = _attach_target_genes(pd.read_csv(path))
    result["condition"] = "EBs"
    return result


@st.cache_data(ttl=3600)
def load_deg_table(perturbation: str) -> Optional[pd.DataFrame]:
    """Load an EBs differential-expression table for one perturbation."""
    if not DEG_TABLES_EBS.exists():
        return None
    path = DEG_TABLES_EBS / f"{perturbation}_deg.csv"
    if path.exists():
        return pd.read_csv(path)
    matches = list(DEG_TABLES_EBS.glob(f"{perturbation}*_deg.csv"))
    return pd.read_csv(matches[0]) if matches else None


@st.cache_data(ttl=3600)
def load_dose_response() -> Optional[pd.DataFrame]:
    """Load EBs dose-response summaries."""
    path = DOSE_RESPONSE_EBS / "EBs_dose_response.csv"
    return _attach_target_genes(pd.read_csv(path)) if path.exists() else None


@st.cache_data(ttl=3600)
def load_pathway_enrichment_for_perturbation(
    perturbation: str,
) -> Optional[pd.DataFrame]:
    """Load pathway-enrichment rows for one EBs perturbation.

    The filtered Parquet export supports predicate pushdown. The chunked CSV
    fallback bounds memory use for older local data bundles.
    """
    parquet_path = DATA_EXTRACTED / "pathway_enrichment_ebs.parquet"
    if parquet_path.exists():
        return pd.read_parquet(
            parquet_path, filters=[("perturbation", "==", perturbation)]
        )

    path = PATHWAY_ENRICHMENT_EBS / "EBs_enrichment_results.csv"
    if not path.exists():
        return None
    matches = []
    for chunk in pd.read_csv(path, chunksize=100_000):
        selected = chunk[chunk["perturbation"] == perturbation]
        if not selected.empty:
            matches.append(selected)
    return pd.concat(matches, ignore_index=True) if matches else None


@st.cache_data(ttl=3600)
def load_transcriptome_edist() -> Optional[pd.DataFrame]:
    """Load EBs transcriptome E-distance summaries."""
    path = TRANSCRIPTOME_EDIST_EBS / "EBs_tedist_merged.csv"
    return _attach_target_genes(pd.read_csv(path)) if path.exists() else None


@st.cache_data(ttl=3600)
def load_tf_clusters() -> Optional[pd.DataFrame]:
    """Load EBs transcriptomic-signature cluster assignments."""
    path = TF_SIMILARITY_EBS / "EBs_tf_clusters.csv"
    return pd.read_csv(path) if path.exists() else None


@st.cache_data(ttl=3600)
def load_lineage_de(
    perturbation: str, lineage: str
) -> Optional[pd.DataFrame]:
    """Load an EBs lineage-specific differential-expression table."""
    path = LINEAGE_DE_EBS / f"{perturbation}_{lineage}_DE.csv"
    return pd.read_csv(path) if path.exists() else None


def get_available_lineage_de(perturbation: str) -> list[str]:
    """Return lineages with an exported DE table for a perturbation."""
    if not LINEAGE_DE_EBS.exists():
        return []
    lineages = []
    for path in LINEAGE_DE_EBS.glob(f"{perturbation}_*_DE.csv"):
        suffix = path.stem[len(perturbation) + 1 :]
        if suffix.endswith("_DE"):
            lineages.append(suffix[:-3])
    return sorted(lineages)


@st.cache_data(ttl=3600)
def load_timecourse_expression() -> pd.DataFrame:
    """Load gene-level baseline timecourse expression and detection summaries."""
    path = DATA_EXTRACTED / "timecourse_expression.parquet"
    if not path.exists():
        raise FileNotFoundError(
            "Timecourse expression is not extracted; run "
            "scripts/extract_timecourse_expression.py"
        )
    return pd.read_parquet(path)


@st.cache_data(ttl=3600)
def get_gene_list() -> list[str]:
    """Return canonical symbols that have a usable EBs phenotype profile."""
    lineage = load_lineage_analysis()
    return sorted(lineage["gene"].dropna().astype(str).unique().tolist())


@st.cache_data(ttl=3600)
def get_perturbations_for_gene(gene: str) -> dict:
    """Return all EBs perturbations mapped to a canonical target symbol."""
    lineage = load_lineage_analysis()
    perturbations = (
        lineage.loc[lineage["gene"] == gene, "perturbation"].drop_duplicates().tolist()
    )
    return {
        "perturbations": perturbations,
        "conditions": {value: ["EBs"] for value in perturbations},
    }


@st.cache_data(ttl=3600)
def get_summary_stats() -> dict:
    """Derive release summary statistics from the current exports."""
    try:
        lineage = load_lineage_analysis()
    except FileNotFoundError as exc:
        return {"error": str(exc)}

    unique_perturbations = lineage.drop_duplicates("perturbation")
    stats = {
        "n_genes": int(lineage["gene"].nunique()),
        "n_perturbations": int(lineage["perturbation"].nunique()),
        "n_cells_ebs": int(unique_perturbations["n_cells"].sum()),
        "n_lineages": int(lineage["lineage"].nunique()),
        "n_cells_timecourse": 0,
        "n_timepoints": 0,
        "n_screen_genes_timecourse": 0,
    }
    try:
        timecourse = load_timecourse_expression()
        stats["n_cells_timecourse"] = int(
            timecourse.drop_duplicates("day")["n_cells"].sum()
        )
        stats["n_timepoints"] = int(timecourse["day"].nunique())
        stats["n_screen_genes_timecourse"] = len(
            set(lineage["gene"].dropna().astype(str))
            & set(timecourse["gene"].dropna().astype(str))
        )
    except FileNotFoundError:
        pass
    return stats


@st.cache_data(ttl=3600)
def get_data_status() -> dict:
    """Return human-readable provenance for the data shown in the portal."""
    core_path = RESULTS_EBS / "lineage_analysis" / "EBs_lineage_analysis_merged.csv"
    status = {
        "scope": "Embryoid bodies (EBs)",
        "screen_updated": None,
        "timecourse_updated": None,
        "timecourse_scale": None,
    }
    if core_path.exists():
        status["screen_updated"] = datetime.fromtimestamp(
            core_path.stat().st_mtime
        ).date().isoformat()

    manifest_path = DATA_EXTRACTED / "timecourse_manifest.json"
    if manifest_path.exists():
        try:
            manifest = json.loads(manifest_path.read_text())
            status["timecourse_updated"] = str(manifest.get("generated_at", ""))[
                :10
            ] or None
            status["timecourse_scale"] = manifest.get("expression_scale")
        except (OSError, ValueError):
            pass

    if status["timecourse_updated"] is None:
        timecourse_path = DATA_EXTRACTED / "timecourse_expression.parquet"
        if timecourse_path.exists():
            status["timecourse_updated"] = datetime.fromtimestamp(
                timecourse_path.stat().st_mtime
            ).date().isoformat()
    return status
