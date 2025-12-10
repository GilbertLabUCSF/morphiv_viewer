#!/bin/bash
# MORPHIC Portal - Setup and Run Script
set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

ENV_DIR="./morphic_website_env"
PORT="${1:-8501}"

echo "=== MORPHIC TF Perturbation Portal ==="
echo

# Check environment
if [ ! -d "$ENV_DIR" ]; then
    echo "ERROR: Virtual environment not found at $ENV_DIR"
    echo "Create it with: python -m venv morphic_website_env && pip install -r requirements.txt"
    exit 1
fi

# Activate environment
source "$ENV_DIR/bin/activate"

# Validate data sources
echo "[1/3] Validating data sources..."
python scripts/validate_data_sources.py
echo

# Extract timecourse data if needed
TIMECOURSE_FILE="data_extracted/timecourse_expression.parquet"
if [ ! -f "$TIMECOURSE_FILE" ]; then
    echo "[2/3] Extracting timecourse expression data..."
    python scripts/extract_timecourse_expression.py
    echo
else
    echo "[2/3] Timecourse data already extracted, skipping"
    echo
fi

# Launch server
echo "[3/3] Starting Streamlit on port $PORT..."
echo "      Access at: http://localhost:$PORT"
echo
exec "$ENV_DIR/bin/streamlit" run app.py --server.headless true --server.port "$PORT"
