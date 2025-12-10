#!/bin/bash
# Upload all data files to GCP instance for full portal functionality
# Run from morphic_website directory
# Total size: ~240 MB

set -e

SOURCE_DATA="/large_storage/gilbertlab/ashir/morphic_pub_refactor/tf_perturbseq"
GCP_INSTANCE="instance-20251119-174124"
GCP_ZONE="us-central1-c"
GCP_PROJECT="ashir-borah-project"
REMOTE_DATA="/opt/morphic/data/tf_perturbseq"

echo "=== MORPHIC Portal Data Upload ==="
echo "Total data: ~240 MB"
echo ""

# Create all required directories
echo "[1/7] Creating remote directories..."
gcloud compute ssh --zone "$GCP_ZONE" "$GCP_INSTANCE" --project "$GCP_PROJECT" -- "
sudo mkdir -p $REMOTE_DATA/results/EBs/{lineage_analysis,knockdown_efficiency,compositional}
sudo mkdir -p $REMOTE_DATA/results/iPSC/{lineage_analysis,knockdown_efficiency,compositional}
sudo mkdir -p $REMOTE_DATA/results_dec1/EBs/viability
sudo mkdir -p $REMOTE_DATA/results_dec1/iPSC/viability
sudo mkdir -p $REMOTE_DATA/figures
sudo mkdir -p /opt/morphic/website/data_extracted
sudo chown -R \$USER /opt/morphic
"

# Upload lineage analysis (~20 MB)
echo "[2/7] Uploading lineage analysis (~20 MB)..."
gcloud compute scp --zone "$GCP_ZONE" --project "$GCP_PROJECT" \
  "$SOURCE_DATA/results/EBs/lineage_analysis/EBs_lineage_analysis_merged.csv" \
  "$GCP_INSTANCE:$REMOTE_DATA/results/EBs/lineage_analysis/"
gcloud compute scp --zone "$GCP_ZONE" --project "$GCP_PROJECT" \
  "$SOURCE_DATA/results/iPSC/lineage_analysis/iPSC_lineage_analysis_merged.csv" \
  "$GCP_INSTANCE:$REMOTE_DATA/results/iPSC/lineage_analysis/"

# Upload knockdown efficiency (~700 KB)
echo "[3/7] Uploading knockdown efficiency..."
gcloud compute scp --zone "$GCP_ZONE" --project "$GCP_PROJECT" \
  "$SOURCE_DATA/results/EBs/knockdown_efficiency/knockdown_efficiency_all_genes.csv" \
  "$GCP_INSTANCE:$REMOTE_DATA/results/EBs/knockdown_efficiency/"
gcloud compute scp --zone "$GCP_ZONE" --project "$GCP_PROJECT" \
  "$SOURCE_DATA/results/iPSC/knockdown_efficiency/knockdown_efficiency_all_genes.csv" \
  "$GCP_INSTANCE:$REMOTE_DATA/results/iPSC/knockdown_efficiency/"

# Upload compositional data (~200 KB)
echo "[4/7] Uploading compositional data..."
gcloud compute scp --zone "$GCP_ZONE" --project "$GCP_PROJECT" \
  "$SOURCE_DATA/results/EBs/compositional/EBs_significant_hits.csv" \
  "$GCP_INSTANCE:$REMOTE_DATA/results/EBs/compositional/"
gcloud compute scp --zone "$GCP_ZONE" --project "$GCP_PROJECT" \
  "$SOURCE_DATA/results/iPSC/compositional/iPSC_significant_hits.csv" \
  "$GCP_INSTANCE:$REMOTE_DATA/results/iPSC/compositional/"

# Upload resolved targets and viability (~600 KB)
echo "[5/7] Uploading targets and viability..."
gcloud compute scp --zone "$GCP_ZONE" --project "$GCP_PROJECT" \
  "$SOURCE_DATA/results/EBs/resolved_targets.csv" \
  "$GCP_INSTANCE:$REMOTE_DATA/results/EBs/"
gcloud compute scp --zone "$GCP_ZONE" --project "$GCP_PROJECT" \
  "$SOURCE_DATA/results/iPSC/resolved_targets.csv" \
  "$GCP_INSTANCE:$REMOTE_DATA/results/iPSC/"
gcloud compute scp --zone "$GCP_ZONE" --project "$GCP_PROJECT" \
  "$SOURCE_DATA/results_dec1/EBs/viability/viability_scores_gene_level.csv" \
  "$GCP_INSTANCE:$REMOTE_DATA/results_dec1/EBs/viability/"
gcloud compute scp --zone "$GCP_ZONE" --project "$GCP_PROJECT" \
  "$SOURCE_DATA/results_dec1/iPSC/viability/viability_scores_gene_level.csv" \
  "$GCP_INSTANCE:$REMOTE_DATA/results_dec1/iPSC/viability/"

# Upload figures (~204 MB) - UMAP highlights, dotplots, etc.
echo "[6/7] Uploading figures (~204 MB)..."
gcloud compute scp --zone "$GCP_ZONE" --project "$GCP_PROJECT" --recurse \
  "$SOURCE_DATA/figures/" \
  "$GCP_INSTANCE:$REMOTE_DATA/figures/"

# Upload pre-extracted timecourse data (~2 MB)
echo "[7/7] Uploading timecourse data..."
gcloud compute scp --zone "$GCP_ZONE" --project "$GCP_PROJECT" \
  "data_extracted/timecourse_expression.parquet" \
  "$GCP_INSTANCE:/opt/morphic/website/data_extracted/"

echo ""
echo "=== Upload Complete ==="
echo "All data uploaded to $GCP_INSTANCE"
