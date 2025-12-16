#!/bin/bash
# Quick sync script for iterative development
# Syncs code and extracted data to GCP
# Run from morphic_website directory

set -e

GCLOUD="/home/ashir/google-cloud-sdk/bin/gcloud"
GCP_INSTANCE="instance-20251119-174124"
GCP_ZONE="us-central1-c"
GCP_PROJECT="ashir-borah-project"
REMOTE_DIR="/opt/morphic/website"

echo "=== Syncing MORPHIC Portal to GCP ==="

# Sync Python code
echo "[1/4] Syncing app code..."
$GCLOUD compute scp --zone "$GCP_ZONE" --project "$GCP_PROJECT" \
  app.py \
  "$GCP_INSTANCE:$REMOTE_DIR/"

$GCLOUD compute scp --zone "$GCP_ZONE" --project "$GCP_PROJECT" --recurse \
  pages \
  "$GCP_INSTANCE:$REMOTE_DIR/"

$GCLOUD compute scp --zone "$GCP_ZONE" --project "$GCP_PROJECT" --recurse \
  src \
  "$GCP_INSTANCE:$REMOTE_DIR/"

# Sync config
echo "[2/4] Syncing config..."
$GCLOUD compute scp --zone "$GCP_ZONE" --project "$GCP_PROJECT" \
  .streamlit/config.toml \
  "$GCP_INSTANCE:$REMOTE_DIR/.streamlit/"

# Sync extracted data
echo "[3/4] Syncing extracted data..."
$GCLOUD compute scp --zone "$GCP_ZONE" --project "$GCP_PROJECT" --recurse \
  data_extracted \
  "$GCP_INSTANCE:$REMOTE_DIR/"

# Fix permissions and restart service
echo "[4/4] Fixing permissions and restarting service..."
$GCLOUD compute ssh --zone "$GCP_ZONE" "$GCP_INSTANCE" --project "$GCP_PROJECT" \
  --command="sudo chown www-data:www-data /opt/morphic/website/.streamlit/secrets.toml 2>/dev/null; sudo chmod 600 /opt/morphic/website/.streamlit/secrets.toml 2>/dev/null; sudo systemctl restart morphic"

echo ""
echo "=== Sync Complete ==="
echo "Visit: https://perturb.dev/"
