# MORPHIC Portal GCP Deployment

Deploy the MORPHIC TF Perturbation Screen Portal to GCP Compute Engine.

## Server Details

- **Domain**: https://perturb.dev/
- **cellxgene**: https://cellxgene.perturb.dev/
- **IP**: 34.46.167.158
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

Deploy the nginx configuration:

```bash
# Copy the config
sudo cp deploy/nginx.conf /etc/nginx/sites-available/perturb.dev

# Enable the site
sudo ln -sf /etc/nginx/sites-available/perturb.dev /etc/nginx/sites-enabled/

# Disable default site if it conflicts
sudo rm -f /etc/nginx/sites-enabled/default

# Test and reload
sudo nginx -t
sudo systemctl reload nginx
```

### 7. Setup SSL with Certbot

```bash
# Install certbot if not present
sudo apt update
sudo apt install certbot python3-certbot-nginx -y

# Get certificates for all domains
sudo certbot --nginx -d perturb.dev -d www.perturb.dev -d cellxgene.perturb.dev

# Certbot will automatically:
# - Obtain certificates from Let's Encrypt
# - Modify nginx config to add HTTPS
# - Set up auto-renewal
```

### 8. Configure cellxgene Gateway

Find what port cellxgene gateway is running on:

```bash
# Check listening ports
sudo ss -tlnp | grep -E '(python|cellxgene)'

# If not on port 5005, update nginx config:
sudo sed -i 's/127.0.0.1:5005/127.0.0.1:YOUR_PORT/' /etc/nginx/sites-available/perturb.dev
sudo nginx -t && sudo systemctl reload nginx
```

### 9. DNS Configuration

Ensure these DNS records are set (in GoDaddy or your DNS provider):

| Type | Name | Value |
|------|------|-------|
| A | @ | 34.46.167.158 |
| A | www | 34.46.167.158 |
| A | cellxgene | 34.46.167.158 |

### 10. Test

Visit:
- https://perturb.dev/ (Morphic Portal)
- https://cellxgene.perturb.dev/ (cellxgene Gateway)

## Updating

To deploy updates:

```bash
cd /opt/morphic/website
git pull origin deployment
sudo systemctl restart morphic
```

Or use the sync script from local machine:

```bash
bash scripts/sync_to_gcp.sh
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

## Architecture

```
                    Internet
                        │
                        ▼
              ┌─────────────────┐
              │   34.46.167.158 │
              │   (GCP Instance)│
              └────────┬────────┘
                       │
                       ▼
              ┌─────────────────┐
              │     nginx       │
              │   (port 80/443) │
              └────────┬────────┘
                       │
         ┌─────────────┴─────────────┐
         │                           │
         ▼                           ▼
┌─────────────────┐       ┌─────────────────┐
│   perturb.dev   │       │ cellxgene.      │
│                 │       │ perturb.dev     │
│  Streamlit      │       │  cellxgene      │
│  (port 8501)    │       │  gateway        │
└─────────────────┘       └─────────────────┘
```

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

**SSL certificate issues:**
```bash
sudo certbot renew --dry-run
sudo certbot certificates
```

**Nginx config errors:**
```bash
sudo nginx -t
sudo cat /var/log/nginx/error.log
```
