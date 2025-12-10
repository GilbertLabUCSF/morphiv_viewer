"""
Configuration and constants for the MORPHIC Portal.
"""

import os
from pathlib import Path

# Project paths
PROJECT_ROOT = Path(__file__).parent.parent

# DATA_ROOT can be overridden via environment variable for deployment
# Default: ./data (symlink for local dev)
# GCP: Set MORPHIC_DATA_PATH=/opt/morphic/data
DATA_ROOT = Path(os.environ.get("MORPHIC_DATA_PATH", PROJECT_ROOT / "data"))
DATA_EXTRACTED = PROJECT_ROOT / "data_extracted"

# TF Perturbseq paths
TF_PERTURBSEQ = DATA_ROOT / "tf_perturbseq"
RESULTS_EBS = TF_PERTURBSEQ / "results" / "EBs"
RESULTS_IPSC = TF_PERTURBSEQ / "results" / "iPSC"
VIABILITY_EBS = TF_PERTURBSEQ / "results_dec1" / "EBs" / "viability"
VIABILITY_IPSC = TF_PERTURBSEQ / "results_dec1" / "iPSC" / "viability"
DEG_TABLES = TF_PERTURBSEQ / "results_new_GS" / "EBs" / "differential_expression" / "deg_tables"

# Pre-generated figures from Snakemake pipeline
FIGURES_ROOT = TF_PERTURBSEQ / "figures"
FIGURES_DEC1 = TF_PERTURBSEQ / "figures_dec1"

# Timecourse paths
TIMECOURSE = DATA_ROOT / "single_cell_timecourse"

# Lineages tracked in the screen
LINEAGES = [
    "Amnion",
    "Epiblast",
    "Formative_Epiblast",
    "Neural_Ectoderm",
    "Non_neural_Ectoderm",
    "Trophoblast_Like",
]

# Color schemes
LINEAGE_COLORS = {
    "Amnion": "#e74c3c",
    "Epiblast": "#3498db",
    "Formative_Epiblast": "#9b59b6",
    "Neural_Ectoderm": "#2ecc71",
    "Non_neural_Ectoderm": "#f39c12",
    "Trophoblast_Like": "#1abc9c",
}

CONDITION_COLORS = {
    "iPSC": "#3498db",
    "EBs": "#e74c3c",
}

CONDITIONS = ["iPSC", "EBs"]

# Knockdown categories
KNOCKDOWN_CATEGORIES = ["Excellent", "Good", "Moderate", "Poor"]
KNOCKDOWN_COLORS = {
    "Excellent": "#27ae60",
    "Good": "#2ecc71",
    "Moderate": "#f39c12",
    "Poor": "#e74c3c",
}

# Effect size thresholds
EFFECT_SIZES = ["negligible", "small", "medium", "large"]

# External database URL templates
GENECARDS_URL = "https://www.genecards.org/cgi-bin/carddisp.pl?gene={gene}"
DEPMAP_URL = "https://depmap.org/portal/gene/{gene}?tab=overview"
NCBI_GENE_URL = "https://www.ncbi.nlm.nih.gov/gene/?term={gene}"


def get_genecards_url(gene: str) -> str:
    """Generate GeneCards URL for a gene."""
    return GENECARDS_URL.format(gene=gene)


def get_depmap_url(gene: str) -> str:
    """Generate DepMap URL for a gene."""
    return DEPMAP_URL.format(gene=gene)


def get_ncbi_gene_url(gene: str) -> str:
    """Generate NCBI Gene URL for a gene."""
    return NCBI_GENE_URL.format(gene=gene)


def get_umap_highlight_path(perturbation: str, condition: str, fmt: str = "png") -> Path:
    """Get path to pre-generated UMAP highlight figure."""
    return FIGURES_ROOT / condition / "umap" / "highlights" / f"{perturbation}_highlight.{fmt}"


def get_marker_dotplot_path(perturbation: str, condition: str, fmt: str = "png") -> Path:
    """Get path to pre-generated marker dotplot figure."""
    return FIGURES_ROOT / "perturbation_validation" / condition / f"{perturbation}_markers_dotplot.{fmt}"


def get_anchor_dotplot_path(condition: str, variant: str = "standardized", fmt: str = "png") -> Path:
    """Get path to pre-generated anchor dotplot figure.

    Args:
        condition: iPSC or EBs
        variant: standardized, positive, or negative
        fmt: png or svg
    """
    return FIGURES_ROOT / "anchor_validation" / f"{condition}_anchor_dotplot_{variant}.{fmt}"
