#!/bin/bash
# Upload all data files to GCP instance for full portal functionality
# Run from morphic_website directory
# Total size: ~16 GB (includes full DEG tables, lineage DE, and CellxGene h5ad files)

set -e

# GCP SDK path
GCLOUD="/home/ashir/google-cloud-sdk/bin/gcloud"

SOURCE_DATA="/large_storage/gilbertlab/ashir/morphic_pub_refactor/tf_perturbseq"
GCP_INSTANCE="instance-20251119-174124"
GCP_ZONE="us-central1-c"
GCP_PROJECT="ashir-borah-project"
REMOTE_DATA="/opt/morphic/data/tf_perturbseq"

echo "=== MORPHIC Portal Data Upload ==="
echo "Total data: ~16 GB"
echo ""

# Create all required directories
echo "[1/12] Creating remote directories..."
$GCLOUD compute ssh --zone "$GCP_ZONE" "$GCP_INSTANCE" --project "$GCP_PROJECT" -- "
sudo mkdir -p $REMOTE_DATA/results/EBs/{lineage_analysis,knockdown_efficiency,compositional,viability,differential_expression/deg_tables,dose_response,pathway_enrichment,tf_similarity,transcriptome_edist,lineage_de}
sudo mkdir -p $REMOTE_DATA/results/iPSC/{lineage_analysis,knockdown_efficiency,compositional,viability,differential_expression/deg_tables,dose_response,pathway_enrichment,tf_similarity,transcriptome_edist,lineage_de}
sudo mkdir -p $REMOTE_DATA/figures/{EBs,iPSC}/{lineage_analysis/spider,dose_response,pathway_enrichment,tf_similarity,transcriptome_edist}
sudo mkdir -p /opt/morphic/website/data_extracted
sudo mkdir -p /opt/morphic/cellxgene
sudo chown -R \$USER /opt/morphic
"

# Upload lineage analysis (~20 MB)
echo "[2/12] Uploading lineage analysis (~20 MB)..."
$GCLOUD compute scp --zone "$GCP_ZONE" --project "$GCP_PROJECT" \
  "$SOURCE_DATA/results/EBs/lineage_analysis/EBs_lineage_analysis_merged.csv" \
  "$GCP_INSTANCE:$REMOTE_DATA/results/EBs/lineage_analysis/"
$GCLOUD compute scp --zone "$GCP_ZONE" --project "$GCP_PROJECT" \
  "$SOURCE_DATA/results/iPSC/lineage_analysis/iPSC_lineage_analysis_merged.csv" \
  "$GCP_INSTANCE:$REMOTE_DATA/results/iPSC/lineage_analysis/"

# Upload knockdown efficiency (~700 KB)
echo "[3/12] Uploading knockdown efficiency..."
$GCLOUD compute scp --zone "$GCP_ZONE" --project "$GCP_PROJECT" \
  "$SOURCE_DATA/results/EBs/knockdown_efficiency/knockdown_efficiency_all_genes.csv" \
  "$GCP_INSTANCE:$REMOTE_DATA/results/EBs/knockdown_efficiency/"
$GCLOUD compute scp --zone "$GCP_ZONE" --project "$GCP_PROJECT" \
  "$SOURCE_DATA/results/iPSC/knockdown_efficiency/knockdown_efficiency_all_genes.csv" \
  "$GCP_INSTANCE:$REMOTE_DATA/results/iPSC/knockdown_efficiency/"

# Upload compositional data (~200 KB)
echo "[4/12] Uploading compositional data..."
$GCLOUD compute scp --zone "$GCP_ZONE" --project "$GCP_PROJECT" \
  "$SOURCE_DATA/results/EBs/compositional/EBs_significant_hits.csv" \
  "$GCP_INSTANCE:$REMOTE_DATA/results/EBs/compositional/"
$GCLOUD compute scp --zone "$GCP_ZONE" --project "$GCP_PROJECT" \
  "$SOURCE_DATA/results/iPSC/compositional/iPSC_significant_hits.csv" \
  "$GCP_INSTANCE:$REMOTE_DATA/results/iPSC/compositional/"

# Upload resolved targets and viability (~600 KB)
echo "[5/12] Uploading targets and viability..."
for COND in EBs iPSC; do
    $GCLOUD compute scp --zone "$GCP_ZONE" --project "$GCP_PROJECT" \
      "$SOURCE_DATA/results/$COND/resolved_targets.csv" \
      "$GCP_INSTANCE:$REMOTE_DATA/results/$COND/"
    $GCLOUD compute scp --zone "$GCP_ZONE" --project "$GCP_PROJECT" \
      "$SOURCE_DATA/results/$COND/viability/viability_scores_gene_level.csv" \
      "$GCP_INSTANCE:$REMOTE_DATA/results/$COND/viability/"
done

# Upload figures (~204 MB)
echo "[6/12] Uploading figures (~204 MB)..."
$GCLOUD compute scp --zone "$GCP_ZONE" --project "$GCP_PROJECT" --recurse \
  "$SOURCE_DATA/figures/EBs" \
  "$SOURCE_DATA/figures/iPSC" \
  "$SOURCE_DATA/figures/anchor_validation" \
  "$SOURCE_DATA/figures/knockdown_efficiency" \
  "$SOURCE_DATA/figures/perturbation_validation" \
  "$SOURCE_DATA/figures/supplementary" \
  "$GCP_INSTANCE:$REMOTE_DATA/figures/"

# Upload pre-extracted timecourse data (~2 MB)
echo "[7/12] Uploading timecourse data..."
$GCLOUD compute scp --zone "$GCP_ZONE" --project "$GCP_PROJECT" \
  "data_extracted/timecourse_expression.parquet" \
  "$GCP_INSTANCE:/opt/morphic/website/data_extracted/"

# Upload ALL DEG tables (~6.4 GB per condition)
echo "[8/12] Uploading DEG tables (EBs + iPSC, ~12.8 GB total)..."
for COND in EBs iPSC; do
    DEG_DIR="$SOURCE_DATA/results/$COND/differential_expression/deg_tables"
    if [ -d "$DEG_DIR" ]; then
        $GCLOUD compute scp --zone "$GCP_ZONE" --project "$GCP_PROJECT" --recurse \
          "$DEG_DIR" \
          "$GCP_INSTANCE:$REMOTE_DATA/results/$COND/differential_expression/"
    fi
done

# Upload new small CSVs (~175 MB total)
echo "[9/12] Uploading new analysis data (~175 MB)..."
for COND in EBs iPSC; do
    # Dose response (~68 KB)
    DR_FILE="$SOURCE_DATA/results/$COND/dose_response/${COND}_dose_response.csv"
    if [ -f "$DR_FILE" ]; then
        $GCLOUD compute scp --zone "$GCP_ZONE" --project "$GCP_PROJECT" \
          "$DR_FILE" \
          "$GCP_INSTANCE:$REMOTE_DATA/results/$COND/dose_response/"
    fi

    # Pathway enrichment (~125 MB)
    PE_FILE="$SOURCE_DATA/results/$COND/pathway_enrichment/${COND}_enrichment_results.csv"
    if [ -f "$PE_FILE" ]; then
        $GCLOUD compute scp --zone "$GCP_ZONE" --project "$GCP_PROJECT" \
          "$PE_FILE" \
          "$GCP_INSTANCE:$REMOTE_DATA/results/$COND/pathway_enrichment/"
    fi

    # TF similarity (~48 MB)
    for TF_FILE in "${COND}_tf_clusters.csv" "${COND}_tf_similarity_matrix.csv"; do
        TF_PATH="$SOURCE_DATA/results/$COND/tf_similarity/$TF_FILE"
        if [ -f "$TF_PATH" ]; then
            $GCLOUD compute scp --zone "$GCP_ZONE" --project "$GCP_PROJECT" \
              "$TF_PATH" \
              "$GCP_INSTANCE:$REMOTE_DATA/results/$COND/tf_similarity/"
        fi
    done

    # Transcriptome E-distance (~452 KB)
    ED_FILE="$SOURCE_DATA/results/$COND/transcriptome_edist/${COND}_tedist_merged.csv"
    if [ -f "$ED_FILE" ]; then
        $GCLOUD compute scp --zone "$GCP_ZONE" --project "$GCP_PROJECT" \
          "$ED_FILE" \
          "$GCP_INSTANCE:$REMOTE_DATA/results/$COND/transcriptome_edist/"
    fi
done

# Upload ALL lineage DE files (~8.1 GB per condition)
echo "[10/12] Uploading lineage DE files (EBs + iPSC)..."
for COND in EBs iPSC; do
    LDE_DIR="$SOURCE_DATA/results/$COND/lineage_de"
    if [ -d "$LDE_DIR" ]; then
        $GCLOUD compute scp --zone "$GCP_ZONE" --project "$GCP_PROJECT" --recurse \
          "$LDE_DIR" \
          "$GCP_INSTANCE:$REMOTE_DATA/results/$COND/"
    fi
done

# Upload CellxGene h5ad files (26GB EBs + 15GB iPSC)
echo "[11/12] Uploading CellxGene h5ad files (~41 GB)..."
for COND in EBs iPSC; do
    H5AD_FILE="$SOURCE_DATA/results/$COND/cell_type_scored_scvi_compact.h5ad"
    if [ -f "$H5AD_FILE" ]; then
        echo "  Uploading ${COND} h5ad ($(du -h "$H5AD_FILE" | cut -f1))..."
        $GCLOUD compute scp --zone "$GCP_ZONE" --project "$GCP_PROJECT" \
          "$H5AD_FILE" \
          "$GCP_INSTANCE:/opt/morphic/cellxgene/${COND}_cell_type_scored_scvi_compact.h5ad"
    fi
done

# Configure CellxGene Gateway to serve new files
echo "[12/12] Configuring CellxGene..."
$GCLOUD compute ssh --zone "$GCP_ZONE" "$GCP_INSTANCE" --project "$GCP_PROJECT" -- "
# Remove old h5ad files if they exist
rm -f /opt/morphic/cellxgene/nov20_*.h5ad 2>/dev/null || true

# Restart cellxgene gateway if running
if systemctl is-active --quiet cellxgene-gateway; then
    sudo systemctl restart cellxgene-gateway
    echo 'CellxGene gateway restarted'
fi
"

echo ""
echo "=== Upload Complete ==="
echo "All data uploaded to $GCP_INSTANCE"
echo ""
echo "Don't forget to update .streamlit/secrets.toml with the new CellxGene URLs:"
echo '  ipsc_url = "https://cellxgene.perturb.dev/view/iPSC_cell_type_scored_scvi_compact.h5ad/"'
echo '  ebs_url = "https://cellxgene.perturb.dev/view/EBs_cell_type_scored_scvi_compact.h5ad/"'
