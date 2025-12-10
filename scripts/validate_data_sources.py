#!/usr/bin/env python3
"""
Validate all required data sources exist before running the MORPHIC portal.

Run this script before starting the Streamlit app to ensure all data is accessible.
"""

from pathlib import Path
import sys

# Base paths
SCRIPT_DIR = Path(__file__).parent
PROJECT_ROOT = SCRIPT_DIR.parent
DATA_SYMLINK = PROJECT_ROOT / "data"
DATA_EXTRACTED = PROJECT_ROOT / "data_extracted"

# Expected source paths (relative to data symlink)
REQUIRED_FILES = {
    # Lineage analysis
    "tf_perturbseq/results/EBs/lineage_analysis/EBs_lineage_analysis_merged.csv": "EBs lineage analysis",
    "tf_perturbseq/results/iPSC/lineage_analysis/iPSC_lineage_analysis_merged.csv": "iPSC lineage analysis",

    # Knockdown efficiency
    "tf_perturbseq/results/EBs/knockdown_efficiency/knockdown_efficiency_all_genes.csv": "EBs knockdown efficiency",
    "tf_perturbseq/results/iPSC/knockdown_efficiency/knockdown_efficiency_all_genes.csv": "iPSC knockdown efficiency",

    # Viability
    "tf_perturbseq/results_dec1/EBs/viability/viability_scores_gene_level.csv": "EBs viability scores",
    "tf_perturbseq/results_dec1/iPSC/viability/viability_scores_gene_level.csv": "iPSC viability scores",

    # Compositional analysis
    "tf_perturbseq/results/EBs/compositional/EBs_significant_hits.csv": "EBs compositional hits",
    "tf_perturbseq/results/iPSC/compositional/iPSC_significant_hits.csv": "iPSC compositional hits",

    # TF similarity
    "tf_perturbseq/results/EBs/tf_similarity/EBs_tf_clusters.csv": "EBs TF clusters",

    # Gene list
    "tf_perturbseq/results/EBs/resolved_targets.csv": "EBs resolved targets",
    "tf_perturbseq/results/iPSC/resolved_targets.csv": "iPSC resolved targets",

    # Timecourse h5ad (source for extraction)
    "single_cell_timecourse/results_latest/processed/timecourse_scvi.h5ad": "Timecourse h5ad",
}

REQUIRED_DIRECTORIES = {
    "tf_perturbseq/results_new_GS/EBs/differential_expression/deg_tables": "DEG tables directory",
}

EXTRACTED_FILES = {
    "timecourse_expression.parquet": "Timecourse per-gene expression",
}


def check_symlink():
    """Check that data symlink exists and points to correct location."""
    if not DATA_SYMLINK.exists():
        return False, "Data symlink does not exist"

    if not DATA_SYMLINK.is_symlink():
        return False, "data/ is not a symlink"

    target = DATA_SYMLINK.resolve()
    if not target.exists():
        return False, f"Symlink target does not exist: {target}"

    if "morphic_pub_refactor" not in str(target):
        return False, f"Symlink points to unexpected location: {target}"

    return True, f"Symlink OK: data -> {target}"


def check_required_files():
    """Check all required source files exist."""
    missing = []
    found = []

    for rel_path, description in REQUIRED_FILES.items():
        full_path = DATA_SYMLINK / rel_path
        if full_path.exists():
            found.append((rel_path, description))
        else:
            missing.append((rel_path, description))

    return found, missing


def check_required_directories():
    """Check required directories exist and have content."""
    issues = []
    ok = []

    for rel_path, description in REQUIRED_DIRECTORIES.items():
        full_path = DATA_SYMLINK / rel_path
        if not full_path.exists():
            issues.append((rel_path, description, "Directory does not exist"))
        elif not full_path.is_dir():
            issues.append((rel_path, description, "Path is not a directory"))
        else:
            file_count = len(list(full_path.glob("*.csv")))
            if file_count == 0:
                issues.append((rel_path, description, "Directory is empty"))
            else:
                ok.append((rel_path, description, f"{file_count} files"))

    return ok, issues


def check_extracted_files():
    """Check extracted data files exist."""
    missing = []
    found = []

    for filename, description in EXTRACTED_FILES.items():
        full_path = DATA_EXTRACTED / filename
        if full_path.exists():
            found.append((filename, description))
        else:
            missing.append((filename, description))

    return found, missing


def main():
    print("=" * 60)
    print("MORPHIC Portal Data Validation")
    print("=" * 60)
    print()

    all_ok = True

    # Check symlink
    print("[1/4] Checking data symlink...")
    symlink_ok, symlink_msg = check_symlink()
    if symlink_ok:
        print(f"  ✓ {symlink_msg}")
    else:
        print(f"  ✗ {symlink_msg}")
        print()
        print("  To fix, run:")
        print(f"    cd {PROJECT_ROOT}")
        print("    ln -s ../morphic_pub_refactor data")
        all_ok = False
        print()
        print("Cannot continue without data symlink. Exiting.")
        sys.exit(1)
    print()

    # Check required files
    print("[2/4] Checking required source files...")
    found_files, missing_files = check_required_files()
    for rel_path, description in found_files:
        print(f"  ✓ {description}")
    for rel_path, description in missing_files:
        print(f"  ✗ MISSING: {description}")
        print(f"      Expected: {DATA_SYMLINK / rel_path}")
        all_ok = False
    print()

    # Check required directories
    print("[3/4] Checking required directories...")
    ok_dirs, issue_dirs = check_required_directories()
    for rel_path, description, info in ok_dirs:
        print(f"  ✓ {description} ({info})")
    for rel_path, description, issue in issue_dirs:
        print(f"  ✗ {description}: {issue}")
        print(f"      Expected: {DATA_SYMLINK / rel_path}")
        all_ok = False
    print()

    # Check extracted files
    print("[4/4] Checking extracted data files...")
    found_extracted, missing_extracted = check_extracted_files()
    for filename, description in found_extracted:
        print(f"  ✓ {description}")
    for filename, description in missing_extracted:
        print(f"  ⚠ NOT YET EXTRACTED: {description}")
        print(f"      Run: python scripts/extract_timecourse_expression.py")
    print()

    # Summary
    print("=" * 60)
    if all_ok and not missing_extracted:
        print("✓ All data sources validated successfully!")
        print("  You can now run: streamlit run app.py")
    elif all_ok:
        print("⚠ Source data OK, but extracted data missing.")
        print("  Run extraction scripts before starting the app:")
        print("    python scripts/extract_timecourse_expression.py")
    else:
        print("✗ Some required data sources are missing.")
        print("  Please check the paths above and fix any issues.")
        sys.exit(1)
    print("=" * 60)


if __name__ == "__main__":
    main()
