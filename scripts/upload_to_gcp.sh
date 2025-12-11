#!/bin/bash
# Upload all data files to GCP instance for full portal functionality
# Run from morphic_website directory
# Total size: ~240 MB

set -e

# GCP SDK path
GCLOUD="/home/ashir/google-cloud-sdk/bin/gcloud"

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
$GCLOUD compute ssh --zone "$GCP_ZONE" "$GCP_INSTANCE" --project "$GCP_PROJECT" -- "
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
$GCLOUD compute scp --zone "$GCP_ZONE" --project "$GCP_PROJECT" \
  "$SOURCE_DATA/results/EBs/lineage_analysis/EBs_lineage_analysis_merged.csv" \
  "$GCP_INSTANCE:$REMOTE_DATA/results/EBs/lineage_analysis/"
$GCLOUD compute scp --zone "$GCP_ZONE" --project "$GCP_PROJECT" \
  "$SOURCE_DATA/results/iPSC/lineage_analysis/iPSC_lineage_analysis_merged.csv" \
  "$GCP_INSTANCE:$REMOTE_DATA/results/iPSC/lineage_analysis/"

# Upload knockdown efficiency (~700 KB)
echo "[3/7] Uploading knockdown efficiency..."
$GCLOUD compute scp --zone "$GCP_ZONE" --project "$GCP_PROJECT" \
  "$SOURCE_DATA/results/EBs/knockdown_efficiency/knockdown_efficiency_all_genes.csv" \
  "$GCP_INSTANCE:$REMOTE_DATA/results/EBs/knockdown_efficiency/"
$GCLOUD compute scp --zone "$GCP_ZONE" --project "$GCP_PROJECT" \
  "$SOURCE_DATA/results/iPSC/knockdown_efficiency/knockdown_efficiency_all_genes.csv" \
  "$GCP_INSTANCE:$REMOTE_DATA/results/iPSC/knockdown_efficiency/"

# Upload compositional data (~200 KB)
echo "[4/7] Uploading compositional data..."
$GCLOUD compute scp --zone "$GCP_ZONE" --project "$GCP_PROJECT" \
  "$SOURCE_DATA/results/EBs/compositional/EBs_significant_hits.csv" \
  "$GCP_INSTANCE:$REMOTE_DATA/results/EBs/compositional/"
$GCLOUD compute scp --zone "$GCP_ZONE" --project "$GCP_PROJECT" \
  "$SOURCE_DATA/results/iPSC/compositional/iPSC_significant_hits.csv" \
  "$GCP_INSTANCE:$REMOTE_DATA/results/iPSC/compositional/"

# Upload resolved targets and viability (~600 KB)
echo "[5/7] Uploading targets and viability..."
$GCLOUD compute scp --zone "$GCP_ZONE" --project "$GCP_PROJECT" \
  "$SOURCE_DATA/results/EBs/resolved_targets.csv" \
  "$GCP_INSTANCE:$REMOTE_DATA/results/EBs/"
$GCLOUD compute scp --zone "$GCP_ZONE" --project "$GCP_PROJECT" \
  "$SOURCE_DATA/results/iPSC/resolved_targets.csv" \
  "$GCP_INSTANCE:$REMOTE_DATA/results/iPSC/"
$GCLOUD compute scp --zone "$GCP_ZONE" --project "$GCP_PROJECT" \
  "$SOURCE_DATA/results_dec1/EBs/viability/viability_scores_gene_level.csv" \
  "$GCP_INSTANCE:$REMOTE_DATA/results_dec1/EBs/viability/"
$GCLOUD compute scp --zone "$GCP_ZONE" --project "$GCP_PROJECT" \
  "$SOURCE_DATA/results_dec1/iPSC/viability/viability_scores_gene_level.csv" \
  "$GCP_INSTANCE:$REMOTE_DATA/results_dec1/iPSC/viability/"

# Upload figures (~204 MB) - UMAP highlights, dotplots, etc.
# Note: Using rsync-like behavior to avoid nested directories
echo "[6/8] Uploading figures (~204 MB)..."
$GCLOUD compute scp --zone "$GCP_ZONE" --project "$GCP_PROJECT" --recurse \
  "$SOURCE_DATA/figures/EBs" \
  "$SOURCE_DATA/figures/iPSC" \
  "$SOURCE_DATA/figures/anchor_validation" \
  "$SOURCE_DATA/figures/knockdown_efficiency" \
  "$SOURCE_DATA/figures/perturbation_validation" \
  "$SOURCE_DATA/figures/supplementary" \
  "$GCP_INSTANCE:$REMOTE_DATA/figures/"

# Upload pre-extracted timecourse data (~2 MB)
echo "[7/8] Uploading timecourse data..."
$GCLOUD compute scp --zone "$GCP_ZONE" --project "$GCP_PROJECT" \
  "data_extracted/timecourse_expression.parquet" \
  "$GCP_INSTANCE:/opt/morphic/website/data_extracted/"

# Upload DEG tables for perturbations with UMAP highlights (~240 MB)
echo "[8/8] Uploading DEG tables for highlighted perturbations (~240 MB)..."
$GCLOUD compute ssh --zone "$GCP_ZONE" "$GCP_INSTANCE" --project "$GCP_PROJECT" -- \
  "sudo mkdir -p $REMOTE_DATA/results_new_GS/EBs/differential_expression/deg_tables && sudo chown -R \$USER /opt/morphic"

# List of perturbations with UMAP highlights (EBs)
DEG_PERTURBATIONS="ARID2_P1P2 BCL6_P1P2 C1orf85_P2 CRTC3_P1P2 CRX_P1P2 CSDC2_P1P2 CSDE1_P1P2 CTBP2_P1P2 CTNNB1_P1P2 DLX2_P1P2 DLX3_P1P2 DMBX1_P1P2 DMRTA2_P1P2 DMRTB1_P1P2 DNMT1_P1P2 DNMT3A_P1P2 DPPA2_P1P2 EBF1_P1P2 EBF3_P1P2 ELF1_P1P2 ELK3_P1P2 EP300_P1P2 EP400_P1P2 ETV4_P1P2 ETV5_P1P2 ETV6_P1P2 EZH2_P1 EZH2_P2 FOXD3_P1P2 FOXH1_P1P2 FOXO1_P1 FOXO1_P2 GATA4_P1P2 GBX2_P1P2 GLI2_P1P2 GRHL2_P1P2 GRHL3_P1P2 HES1_P1P2 HES5_P1P2 HES7_P1P2 HNF1A_P1P2 HOXA1_P1P2 HOXB1_P1P2 IRF2_P1P2 KDM1A_P1P2 KDM1B_P1P2 LHX5_P1P2 NANOG_P1P2 NKX2-1_P1P2 OTX2_P1P2 PAX6_P1P2 PRDM14_P1P2 SALL2_P1P2 SALL4_P1P2 SETD2_P1P2 SMAD2_P1P2 SOX17_P1P2 SOX2_P1P2 TFAP2A_P1P2 TFAP2C_P1P2 ZIC2_P1P2 ZIC3_P1P2"

for p in $DEG_PERTURBATIONS; do
    deg_file="$SOURCE_DATA/results_new_GS/EBs/differential_expression/deg_tables/${p}_deg.csv"
    if [ -f "$deg_file" ]; then
        $GCLOUD compute scp --zone "$GCP_ZONE" --project "$GCP_PROJECT" \
          "$deg_file" \
          "$GCP_INSTANCE:$REMOTE_DATA/results_new_GS/EBs/differential_expression/deg_tables/"
    fi
done

echo ""
echo "=== Upload Complete ==="
echo "All data uploaded to $GCP_INSTANCE"
