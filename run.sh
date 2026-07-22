#!/bin/bash
# MORPHIC Portal - Setup and Run Script
set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

ENV_DIR="./morphic_website_env"
PORT="${1:-8501}"

echo "=== MORPHIC EBs TF Perturbation Browser ==="
echo

# Check environment
if [ ! -d "$ENV_DIR" ]; then
    echo "ERROR: Virtual environment not found at $ENV_DIR"
    echo "Create it with: python -m venv morphic_website_env && pip install -r requirements.txt"
    exit 1
fi

# Activate environment
source "$ENV_DIR/bin/activate"

# Refresh the derived timecourse data if needed
TIMECOURSE_FILE="data_extracted/timecourse_expression.parquet"
TIMECOURSE_SOURCE="data/single_cell_timecourse/results/processed/timecourse_merged.h5ad"
if [ ! -f "$TIMECOURSE_FILE" ] || { [ -f "$TIMECOURSE_SOURCE" ] && [ "$TIMECOURSE_SOURCE" -nt "$TIMECOURSE_FILE" ]; }; then
    echo "[1/4] Extracting timecourse expression data..."
    python scripts/extract_timecourse_expression.py
    echo
else
    echo "[1/4] Timecourse export is current, skipping"
    echo
fi

# Build the query-optimized pathway table if the upstream CSV changed
PATHWAY_FILE="data_extracted/pathway_enrichment_ebs.parquet"
PATHWAY_SOURCE="data/tf_perturbseq/results/EBs/pathway_enrichment/EBs_enrichment_results.csv"
if [ -f "$PATHWAY_SOURCE" ] && { [ ! -f "$PATHWAY_FILE" ] || [ "$PATHWAY_SOURCE" -nt "$PATHWAY_FILE" ]; }; then
    echo "[2/4] Building pathway lookup..."
    python scripts/extract_pathway_enrichment.py
    echo
else
    echo "[2/4] Pathway lookup is current, skipping"
    echo
fi

# Validate data sources after all derived data is current
echo "[3/4] Validating data sources..."
python scripts/validate_data_sources.py
echo

# Launch server
echo "[4/4] Starting Streamlit on port $PORT..."
echo "      Access at: http://localhost:$PORT"
echo
exec "$ENV_DIR/bin/streamlit" run app.py --server.headless true --server.port "$PORT"
