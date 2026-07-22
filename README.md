# MORPHIC EBs TF Perturbation Browser

An EBs-focused Streamlit browser for connecting CRISPRi transcription-factor
perturbations to lineage phenotypes and unperturbed developmental expression.

The primary profile uses the manuscript pipeline's exported values, thresholds,
and scales as the source of truth. The spider plot and DEG volcano are rendered
interactively from those same exports, with the static pipeline figures retained
as fallbacks. The spider shows the paper's NTC-centered log2 fold-change in mean
KNN lineage probability; complementary Glass's Δ estimates remain available as
a separate analysis. Timecourse expression is independent baseline context and
must not be interpreted as a perturbation effect.

## Quick Start

```bash
./run.sh
```

This refreshes stale derived exports, validates the release data, and launches
the server on port 8501.

To use a different port:
```bash
./run.sh 8701
```

## Manual Setup

If you need to set things up step by step:

```bash
# 1. Activate environment
source morphic_website_env/bin/activate

# 2. Refresh derived data
python scripts/extract_timecourse_expression.py
python scripts/extract_pathway_enrichment.py

# 3. Validate the EBs release bundle
python scripts/validate_data_sources.py

# 4. Launch server
streamlit run app.py --server.headless true --server.port 8501
```

## Data

Data is loaded via the `data` symlink pointing to `../morphic_pub_refactor`. If validation fails, check that the symlink exists and points to the correct location.

- EBs lineage, knockdown, viability, and supporting analysis files come from
  `data/tf_perturbseq/results/EBs/`.
- `data_extracted/timecourse_expression.parquet` is generated from the full
  merged timecourse H5AD and records mean log1p(CP10K) and percent detected.
- `data_extracted/pathway_enrichment_ebs.parquet` is a query-optimized copy of
  the large pathway CSV. Both derived files are rebuilt by `./run.sh` when stale.

## Release checks

```bash
python -m unittest discover -s tests -v
python scripts/validate_data_sources.py
bash -n run.sh
```

## Troubleshooting

- **Port in use**: Change port with `./run.sh 8701` or kill existing process with `fuser -k 8501/tcp`
- **Import errors referencing system Python**: Make sure the venv is activated or use `./morphic_website_env/bin/streamlit` directly
