#!/bin/bash
# Upload all data files to GCP instance using SSH/rsync
# Run from morphic_website directory
# Total size: ~240 MB
#
# Prerequisites:
#   - SSH key configured for GCP instance
#   - Can connect via: ssh 34.46.167.158

set -e

SOURCE_DATA="/large_storage/gilbertlab/ashir/morphic_pub_refactor/tf_perturbseq"
GCP_IP="34.46.167.158"
GCP_USER="${GCP_USER:-$(whoami)}"  # Override with: GCP_USER=username ./script.sh
REMOTE_DATA="/opt/morphic/data/tf_perturbseq"

echo "=== MORPHIC Portal Data Upload (SSH) ==="
echo "Target: $GCP_USER@$GCP_IP"
echo "Total data: ~240 MB"
echo ""

# Test SSH connection
echo "[0/7] Testing SSH connection..."
if ! ssh -o ConnectTimeout=10 "$GCP_USER@$GCP_IP" "echo 'SSH connection successful'"; then
    echo "ERROR: Cannot connect via SSH. Try:"
    echo "  1. ssh $GCP_USER@$GCP_IP  # Test manually"
    echo "  2. GCP_USER=your_username ./scripts/upload_to_gcp_ssh.sh"
    exit 1
fi

# Create all required directories
echo "[1/7] Creating remote directories..."
ssh "$GCP_USER@$GCP_IP" "
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
scp "$SOURCE_DATA/results/EBs/lineage_analysis/EBs_lineage_analysis_merged.csv" \
    "$GCP_USER@$GCP_IP:$REMOTE_DATA/results/EBs/lineage_analysis/"
scp "$SOURCE_DATA/results/iPSC/lineage_analysis/iPSC_lineage_analysis_merged.csv" \
    "$GCP_USER@$GCP_IP:$REMOTE_DATA/results/iPSC/lineage_analysis/"

# Upload knockdown efficiency (~700 KB)
echo "[3/7] Uploading knockdown efficiency..."
scp "$SOURCE_DATA/results/EBs/knockdown_efficiency/knockdown_efficiency_all_genes.csv" \
    "$GCP_USER@$GCP_IP:$REMOTE_DATA/results/EBs/knockdown_efficiency/"
scp "$SOURCE_DATA/results/iPSC/knockdown_efficiency/knockdown_efficiency_all_genes.csv" \
    "$GCP_USER@$GCP_IP:$REMOTE_DATA/results/iPSC/knockdown_efficiency/"

# Upload compositional data (~200 KB)
echo "[4/7] Uploading compositional data..."
scp "$SOURCE_DATA/results/EBs/compositional/EBs_significant_hits.csv" \
    "$GCP_USER@$GCP_IP:$REMOTE_DATA/results/EBs/compositional/"
scp "$SOURCE_DATA/results/iPSC/compositional/iPSC_significant_hits.csv" \
    "$GCP_USER@$GCP_IP:$REMOTE_DATA/results/iPSC/compositional/"

# Upload resolved targets and viability (~600 KB)
echo "[5/7] Uploading targets and viability..."
scp "$SOURCE_DATA/results/EBs/resolved_targets.csv" \
    "$GCP_USER@$GCP_IP:$REMOTE_DATA/results/EBs/"
scp "$SOURCE_DATA/results/iPSC/resolved_targets.csv" \
    "$GCP_USER@$GCP_IP:$REMOTE_DATA/results/iPSC/"
scp "$SOURCE_DATA/results_dec1/EBs/viability/viability_scores_gene_level.csv" \
    "$GCP_USER@$GCP_IP:$REMOTE_DATA/results_dec1/EBs/viability/"
scp "$SOURCE_DATA/results_dec1/iPSC/viability/viability_scores_gene_level.csv" \
    "$GCP_USER@$GCP_IP:$REMOTE_DATA/results_dec1/iPSC/viability/"

# Upload figures (~204 MB) - using rsync for efficiency
echo "[6/7] Uploading figures (~204 MB)..."
rsync -avz --progress "$SOURCE_DATA/figures/" \
    "$GCP_USER@$GCP_IP:$REMOTE_DATA/figures/"

# Upload pre-extracted timecourse data (~2 MB)
echo "[7/7] Uploading timecourse data..."
scp "data_extracted/timecourse_expression.parquet" \
    "$GCP_USER@$GCP_IP:/opt/morphic/website/data_extracted/"

echo ""
echo "=== Upload Complete ==="
echo "All data uploaded to $GCP_IP"
echo ""
echo "Next steps on GCP server:"
echo "  ssh $GCP_USER@$GCP_IP"
echo "  cd /opt/morphic/website"
echo "  source /opt/morphic/venv/bin/activate"
echo "  MORPHIC_DATA_PATH=/opt/morphic/data streamlit run app.py --server.port 8501"
