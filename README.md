# MORPHIC TF Perturbation Portal

Streamlit app for exploring transcription factor perturbation screen data from the MORPHIC consortium.

## Quick Start

```bash
./run.sh
```

This validates data, extracts any needed files, and launches the server on port 8501.

To use a different port:
```bash
./run.sh 8701
```

## Manual Setup

If you need to set things up step by step:

```bash
# 1. Activate environment
source morphic_website_env/bin/activate

# 2. Validate data symlink
python scripts/validate_data_sources.py

# 3. Extract timecourse data (first time only)
python scripts/extract_timecourse_expression.py

# 4. Launch server
streamlit run app.py --server.headless true --server.port 8501
```

## Data

Data is loaded via the `data` symlink pointing to `../morphic_pub_refactor`. If validation fails, check that the symlink exists and points to the correct location.

## Troubleshooting

- **Port in use**: Change port with `./run.sh 8701` or kill existing process with `fuser -k 8501/tcp`
- **Import errors referencing system Python**: Make sure the venv is activated or use `./morphic_website_env/bin/streamlit` directly
