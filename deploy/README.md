# MORPHIC Portal GCP Deployment

Deploy the MORPHIC TF Perturbation Screen Portal to GCP Compute Engine.

## Server Details

- **IP**: 34.46.167.158
- **URL**: http://34.46.167.158/morphic/
- **SSH**: `gcloud compute ssh --zone "us-central1-c" "instance-20251119-174124" --project "ashir-borah-project"`

## Quick Start

### 1. Upload Data (~240 MB)

From your local machine (where the data resides):

```bash
cd /large_storage/gilbertlab/ashir/morphic_website
bash scripts/upload_to_gcp.sh
```

### 2. Setup Server

SSH to the GCP instance:

```bash
gcloud compute ssh --zone "us-central1-c" "instance-20251119-174124" --project "ashir-borah-project"
```

Then run:

```bash
# Clone the repo
sudo mkdir -p /opt/morphic
sudo chown $USER /opt/morphic
cd /opt/morphic
git clone https://github.com/GilbertLabUCSF/morphiv_viewer.git website
cd website
git checkout deployment

# Create virtualenv
python3 -m venv /opt/morphic/venv
source /opt/morphic/venv/bin/activate
pip install -r requirements-lock.txt

# Create data_extracted directory
mkdir -p data_extracted
```

### 3. Copy Secrets

Copy your secrets.toml to the server (contains OpenAI API key):

```bash
# From local machine
gcloud compute scp --zone "us-central1-c" --project "ashir-borah-project" \
  .streamlit/secrets.toml \
  "instance-20251119-174124:/opt/morphic/website/.streamlit/"
```

### 4. Test Manually

On the GCP server:

```bash
cd /opt/morphic/website
source /opt/morphic/venv/bin/activate
export MORPHIC_DATA_PATH=/opt/morphic/data
streamlit run app.py --server.port 8501 --server.headless true
```

Visit http://34.46.167.158:8501 to test (if firewall allows).

### 5. Setup Systemd Service

```bash
# Copy service file
sudo cp deploy/morphic.service /etc/systemd/system/

# Enable and start
sudo systemctl daemon-reload
sudo systemctl enable morphic
sudo systemctl start morphic

# Check status
sudo systemctl status morphic
```

### 6. Configure Nginx

Add the location block from `deploy/nginx.conf` to your Nginx config:

```bash
sudo nano /etc/nginx/sites-available/default
```

Add the contents of `deploy/nginx.conf` inside your server block.

```bash
sudo nginx -t
sudo systemctl reload nginx
```

### 7. Test

Visit: http://34.46.167.158/morphic/

## Updating

To deploy updates:

```bash
cd /opt/morphic/website
git pull origin deployment
sudo systemctl restart morphic
```

## Logs

```bash
# Streamlit logs
sudo journalctl -u morphic -f

# Nginx logs
sudo tail -f /var/log/nginx/access.log
sudo tail -f /var/log/nginx/error.log
```

## Data Files

Total: ~240 MB

| File | Size | Purpose |
|------|------|---------|
| lineage_analysis/*.csv | 19 MB | Main analysis data |
| knockdown_efficiency/*.csv | 650 KB | Knockdown validation |
| compositional/*.csv | 155 KB | Chi-square results |
| viability/*.csv | 520 KB | Fitness scores |
| resolved_targets.csv | 20 KB | Gene list |
| figures/ | 204 MB | UMAP highlights, dotplots |
| timecourse_expression.parquet | 2 MB | Timecourse data |

## Troubleshooting

**Port already in use:**
```bash
sudo lsof -i :8501
sudo kill <PID>
```

**Service won't start:**
```bash
sudo journalctl -u morphic -n 50
```

**Data not loading:**
```bash
ls -la /opt/morphic/data/tf_perturbseq/results/
```
