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

# Viability — now under results/{cond}/viability/
VIABILITY_EBS = RESULTS_EBS / "viability"
VIABILITY_IPSC = RESULTS_IPSC / "viability"

# DEG tables — split per condition
DEG_TABLES_EBS = RESULTS_EBS / "differential_expression" / "deg_tables"
DEG_TABLES_IPSC = RESULTS_IPSC / "differential_expression" / "deg_tables"

# New analysis paths — Dose Response
DOSE_RESPONSE_EBS = RESULTS_EBS / "dose_response"
DOSE_RESPONSE_IPSC = RESULTS_IPSC / "dose_response"

# Pathway Enrichment
PATHWAY_ENRICHMENT_EBS = RESULTS_EBS / "pathway_enrichment"
PATHWAY_ENRICHMENT_IPSC = RESULTS_IPSC / "pathway_enrichment"

# TF Similarity
TF_SIMILARITY_EBS = RESULTS_EBS / "tf_similarity"
TF_SIMILARITY_IPSC = RESULTS_IPSC / "tf_similarity"

# Transcriptome E-distance
TRANSCRIPTOME_EDIST_EBS = RESULTS_EBS / "transcriptome_edist"
TRANSCRIPTOME_EDIST_IPSC = RESULTS_IPSC / "transcriptome_edist"

# Lineage DE
LINEAGE_DE_EBS = RESULTS_EBS / "lineage_de"
LINEAGE_DE_IPSC = RESULTS_IPSC / "lineage_de"

# Pre-generated figures from Snakemake pipeline
FIGURES_ROOT = TF_PERTURBSEQ / "figures"

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
    "Epiblast": "#5a8fa3",
    "Formative_Epiblast": "#9b59b6",
    "Neural_Ectoderm": "#2ecc71",
    "Non_neural_Ectoderm": "#f39c12",
    "Trophoblast_Like": "#1abc9c",
}

CONDITION_COLORS = {
    "iPSC": "#5a9e8f",
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
    """Get path to pre-generated anchor dotplot figure."""
    return FIGURES_ROOT / "anchor_validation" / f"{condition}_anchor_dotplot_{variant}.{fmt}"


def get_spider_plot_path(perturbation: str, condition: str, fmt: str = "png") -> Path:
    """Get the pipeline-generated manuscript spider plot.

    The radial values are NTC-centered log2 fold changes in mean KNN lineage
    probability. Keep this asset as the authoritative browser rendering so the
    portal and paper always use the same plotting implementation.
    """
    return FIGURES_ROOT / condition / "lineage_analysis" / "spider" / f"{perturbation}_knn_prob_shift.{fmt}"


def get_deg_volcano_path(perturbation: str, condition: str, fmt: str = "png") -> Path:
    """Get the pipeline-generated per-perturbation DEG volcano plot."""
    return FIGURES_ROOT / "volcano_plots" / condition / f"{perturbation}_volcano.{fmt}"


def get_marker_tpm_path(perturbation: str, condition: str, fmt: str = "png") -> Path:
    """Get path to a pre-generated marker CPM validation figure.

    The historical ``marker_tpm`` filename is retained for compatibility with
    the analysis output, but the plotted expression unit is CPM.
    """
    return FIGURES_ROOT / condition / "marker_validation" / f"{perturbation}_marker_tpm.{fmt}"


def get_antibody_validation_path(perturbation: str, condition: str, fmt: str = "png") -> Path:
    """Get path to a pre-generated antibody-marker CPM validation figure."""
    return FIGURES_ROOT / condition / "antibody_validation" / f"{perturbation}_antibody_markers.{fmt}"


def get_dose_response_panel_path(condition: str, fmt: str = "png") -> Path:
    """Get path to pre-generated dose response overview panel."""
    return FIGURES_ROOT / condition / "dose_response" / f"dose_response_panel.{fmt}"


def get_pathway_heatmap_path(condition: str, fmt: str = "png") -> Path:
    """Get path to pre-generated pathway enrichment heatmap."""
    return FIGURES_ROOT / condition / "pathway_enrichment" / f"pathway_heatmap.{fmt}"


def get_tf_similarity_path(condition: str, fmt: str = "png") -> Path:
    """Get path to pre-generated TF similarity annotated figure."""
    return FIGURES_ROOT / condition / "tf_similarity" / f"tf_similarity_annotated.{fmt}"


def get_edist_figure_path(condition: str, fmt: str = "png") -> Path:
    """Get the pipeline-generated top-20 transcriptome E-distance figure."""
    return FIGURES_ROOT / condition / "transcriptome_edist" / f"top20_edist_inverted.{fmt}"


def get_edist_histogram_path(condition: str, fmt: str = "png") -> Path:
    """Get the pipeline-generated transcriptome E-distance distribution."""
    return FIGURES_ROOT / condition / "transcriptome_edist" / f"edist_histogram.{fmt}"


def get_viability_figure_path(condition: str, fmt: str = "png") -> Path:
    """Get the pipeline-generated viability volcano plot."""
    return FIGURES_ROOT / "viability" / condition / f"{condition}_volcano_labeled.{fmt}"


def get_trajectory_figure_path(perturbation: str, condition: str, fmt: str = "png") -> Path:
    """Get path to pre-generated trajectory figure."""
    return FIGURES_ROOT / condition / "trajectory" / f"{perturbation}_trajectory.{fmt}"


def _get_results_dir(condition: str) -> Path:
    """Get results directory for a condition."""
    if condition == "EBs":
        return RESULTS_EBS
    elif condition == "iPSC":
        return RESULTS_IPSC
    else:
        raise ValueError(f"Unknown condition: {condition}")
