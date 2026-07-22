# MORPHIC EBs Browser GCP Deployment

Deploy the MORPHIC TF Perturbation Screen Portal to GCP Compute Engine.

## Server Details

- **Domain**: https://perturb.dev/
- **cellxgene**: https://cellxgene.perturb.dev/
- **IP**: 34.46.167.158
- **SSH**: `gcloud compute ssh --zone "us-central1-c" "instance-20251119-174124" --project "ashir-borah-project"`

## Quick Start

### 1. Upload the release

From your local machine (where the data resides):

```bash
cd /large_storage/gilbertlab/ashir/morphic_website
./scripts/upload_to_gcp.sh
```

This is a resumable full rsync. It uploads the application, derived data, EBs
tables, all paper figures, DEG and lineage-DE tables, and the stable EBs, iPSC,
and timecourse CellxGene datasets. Incomplete H5AD files remain resumable in a
dedicated remote staging directory outside the CellxGene browseable tree.

For a quick code, derived-data, and core EBs-table deployment:

```bash
./scripts/sync_to_gcp.sh
```

Both commands validate the remote release and restart the browser. Current EBs
result directories are exact mirrors; superseded remote files are moved into a
timestamped `/opt/morphic/backups/` release backup. A local lock prevents
scheduled and manual transfers from overlapping.

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

### 3. Test Manually

On the GCP server:

```bash
cd /opt/morphic/website
source /opt/morphic/venv/bin/activate
export MORPHIC_DATA_PATH=/opt/morphic/data
streamlit run app.py --server.port 8501 --server.headless true
```

Visit http://34.46.167.158:8501 to test (if firewall allows).

### 4. Setup Systemd Service

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

### 5. Configure Nginx

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

### 6. Setup SSL with Certbot

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

### 7. Configure cellxgene Gateway

Find what port cellxgene gateway is running on:

```bash
# Check listening ports
sudo ss -tlnp | grep -E '(python|cellxgene)'

# If not on port 5000, update nginx config:
sudo sed -i 's/127.0.0.1:5000/127.0.0.1:YOUR_PORT/' /etc/nginx/sites-available/perturb.dev
sudo nginx -t && sudo systemctl reload nginx
```

### 8. DNS Configuration

Ensure these DNS records are set (in GoDaddy or your DNS provider):

| Type | Name | Value |
|------|------|-------|
| A | @ | 34.46.167.158 |
| A | www | 34.46.167.158 |
| A | cellxgene | 34.46.167.158 |

### 9. Test

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

Or deploy code and core EBs data from the local checkout:

```bash
./scripts/sync_to_gcp.sh
```

Use `./scripts/upload_to_gcp.sh` when figures, large analysis tables, or
CellxGene H5AD files have changed. To schedule a nightly full refresh, install
the template in `deploy/morphic-rsync.cron.example` with `crontab -e` after
checking its checkout and log paths.

## Logs

```bash
# Streamlit logs
sudo journalctl -u morphic -f

# Nginx logs
sudo tail -f /var/log/nginx/access.log
sudo tail -f /var/log/nginx/error.log
```

## Data Files

The quick sync includes the query-facing EBs tables below. The full sync also
includes all manuscript figures, DEG and lineage-DE tables, and the three
CellxGene datasets; it is intentionally multi-GB and safely resumable.

| File | Size | Purpose |
|------|------|---------|
| lineage_analysis/*.csv | ~19 MB | Glass's Δ and significance |
| knn_probability_shifts/*.csv | ~3 MB | Interactive spider profiles |
| knockdown_efficiency/*.csv | <1 MB | Knockdown validation |
| compositional/*.csv | ~3 MB | Probability-shift hits |
| viability/*.csv | <1 MB | Fitness scores |
| timecourse_expression.parquet | ~3 MB | Timecourse summary |
| pathway_enrichment_ebs.parquet | ~95 MB | Filterable pathway lookup |

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
