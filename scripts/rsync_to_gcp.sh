#!/bin/bash
# Resumable production sync for the MORPHIC EBs browser.
#
# Usage:
#   scripts/rsync_to_gcp.sh --quick       # code, derived data, core EBs tables
#   scripts/rsync_to_gcp.sh --full        # quick sync + figures, DE, CellxGene
#   scripts/rsync_to_gcp.sh --full --dry-run
#
# Stable CellxGene destination names make repeated runs and cron safe. Rsync
# keeps incomplete H5AD transfers outside the gateway's browseable data tree
# and atomically renames each completed file into place.

set -euo pipefail

MODE="quick"
DRY_RUN=0
for argument in "$@"; do
    case "$argument" in
        --quick) MODE="quick" ;;
        --full) MODE="full" ;;
        --dry-run) DRY_RUN=1 ;;
        *) echo "Unknown argument: $argument"; exit 2 ;;
    esac
done

PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
DATA_ROOT="${MORPHIC_SOURCE_DATA:-$PROJECT_ROOT/data/tf_perturbseq}"
REMOTE_HOST="${MORPHIC_REMOTE_HOST:-ashir@34.46.167.158}"
REMOTE_ROOT="${MORPHIC_REMOTE_ROOT:-/opt/morphic}"
SSH_KEY="${MORPHIC_SSH_KEY:-/home/ashir/.ssh/google_compute_engine}"
LOCK_FILE="${MORPHIC_SYNC_LOCK:-/tmp/morphic-website-rsync.lock}"
CELLXGENE_DATA_DIR="${MORPHIC_CELLXGENE_DATA_DIR:-/home/ashir/cellxgene_data}"
CELLXGENE_PARTIAL_DIR="${MORPHIC_CELLXGENE_PARTIAL_DIR:-/home/ashir/cellxgene_upload_partial}"
RUN_ID="$(date -u +%Y%m%dT%H%M%SZ)"

EB_CELLXGENE="$DATA_ROOT/results/EBs/scanvi_annotated_with_anchors.h5ad"
IPSC_CELLXGENE="$DATA_ROOT/results_20260716/iPSC/scanvi_annotated_with_anchors.h5ad"
TIMECOURSE_CELLXGENE="$PROJECT_ROOT/data/single_cell_timecourse/results/processed/timecourse_scanvi_annotated_with_anchors.h5ad"

SSH_COMMAND="ssh -i $SSH_KEY -o IdentitiesOnly=yes -o BatchMode=yes -o ServerAliveInterval=30 -o ServerAliveCountMax=20"
SSH=(ssh -i "$SSH_KEY" -o IdentitiesOnly=yes -o BatchMode=yes -o ServerAliveInterval=30 -o ServerAliveCountMax=20)
RSYNC_COMMON=(
    rsync
    --archive
    --no-owner
    --no-group
    --human-readable
    --itemize-changes
    --info=progress2,stats2
    --chmod=Du=rwx,Dgo=rx,Fu=rw,Fgo=r
    -e "$SSH_COMMAND"
)
RSYNC=("${RSYNC_COMMON[@]}" --partial --partial-dir=.rsync-partial)
CELLXGENE_RSYNC=(
    "${RSYNC_COMMON[@]}"
    --partial
    --partial-dir="$CELLXGENE_PARTIAL_DIR"
)
if [ "$DRY_RUN" = "1" ]; then
    RSYNC+=(--dry-run)
    CELLXGENE_RSYNC+=(--dry-run)
fi

exec 9>"$LOCK_FILE"
if ! flock -n 9; then
    echo "Another MORPHIC sync is already running ($LOCK_FILE)."
    exit 1
fi

cd "$PROJECT_ROOT"

REQUIRED_PATHS=(
    "$SSH_KEY"
    "$DATA_ROOT/results/EBs"
)
if [ "$MODE" = "full" ]; then
    REQUIRED_PATHS+=(
        "$DATA_ROOT/figures"
        "$DATA_ROOT/results/EBs/differential_expression/deg_tables"
        "$DATA_ROOT/results/EBs/lineage_de"
        "$EB_CELLXGENE"
        "$IPSC_CELLXGENE"
        "$TIMECOURSE_CELLXGENE"
    )
fi

for required in "${REQUIRED_PATHS[@]}"; do
    if [ ! -e "$required" ]; then
        echo "Required path is unavailable: $required"
        exit 1
    fi
done

remote() {
    "${SSH[@]}" "$REMOTE_HOST" "$@"
}

sync_directory() {
    local source="$1"
    local destination="$2"
    if [ ! -d "$source" ]; then
        echo "Skipping unavailable directory: $source"
        return
    fi
    echo
    echo "Syncing $source -> $destination"
    "${RSYNC[@]}" "$source/" "$REMOTE_HOST:$destination/"
}

mirror_directory() {
    local source="$1"
    local destination="$2"
    local backup_label="$3"
    if [ ! -d "$source" ]; then
        echo "Skipping unavailable directory: $source"
        return
    fi
    if [ -z "$(find "$source" -mindepth 1 -print -quit)" ]; then
        echo "Refusing to mirror an empty source directory: $source"
        exit 1
    fi

    local mirror_rsync=("${RSYNC[@]}" --delete-delay)
    if [ "$DRY_RUN" = "0" ]; then
        local backup_dir="$REMOTE_ROOT/backups/$RUN_ID/data/$backup_label"
        remote "mkdir -p $backup_dir"
        mirror_rsync+=(--backup --backup-dir="$backup_dir")
    fi

    echo
    echo "Mirroring $source -> $destination"
    "${mirror_rsync[@]}" "$source/" "$REMOTE_HOST:$destination/"
}

sync_cellxgene_file() {
    local source="$1"
    local destination="$2"
    if [ ! -f "$source" ]; then
        echo "Skipping unavailable file: $source"
        return
    fi
    echo
    echo "Syncing $source -> $destination"
    "${CELLXGENE_RSYNC[@]}" "$source" "$REMOTE_HOST:$destination"
}

echo "=== MORPHIC production rsync ($MODE) ==="
echo "Remote: $REMOTE_HOST · run: $RUN_ID · dry run: $DRY_RUN"

if [ "$DRY_RUN" = "0" ]; then
    remote "mkdir -p \
        $REMOTE_ROOT/backups/$RUN_ID \
        $REMOTE_ROOT/website/data_extracted \
        $REMOTE_ROOT/data/tf_perturbseq/results/EBs \
        $REMOTE_ROOT/data/tf_perturbseq/figures \
        $CELLXGENE_DATA_DIR \
        $CELLXGENE_PARTIAL_DIR"
fi

echo
echo "[1/4] Syncing application code"
CODE_RSYNC=("${RSYNC[@]}" --delete-delay --backup --backup-dir="$REMOTE_ROOT/backups/$RUN_ID")
CODE_RSYNC+=(
    --exclude=.git/
    --exclude=.agents/
    --exclude=.claude/
    --exclude=.codex/
    --exclude=.streamlit/secrets.toml
    --exclude=morphic_website_env/
    --exclude=data
    --exclude=data_extracted/
    --exclude=__pycache__/
    --exclude='*.pyc'
)
"${CODE_RSYNC[@]}" "$PROJECT_ROOT/" "$REMOTE_HOST:$REMOTE_ROOT/website/"
sync_directory "$PROJECT_ROOT/data_extracted" "$REMOTE_ROOT/website/data_extracted"

echo
echo "[2/4] Syncing release-facing EBs tables"
CORE_RESULT_DIRS=(
    lineage_analysis
    knockdown_efficiency
    knn_probability_shifts
    compositional
    viability
    dose_response
    pathway_enrichment
    tf_similarity
    transcriptome_edist
)
for directory in "${CORE_RESULT_DIRS[@]}"; do
    mirror_directory \
        "$DATA_ROOT/results/EBs/$directory" \
        "$REMOTE_ROOT/data/tf_perturbseq/results/EBs/$directory" \
        "core_$directory"
done

if [ "$DRY_RUN" = "0" ]; then
    echo
    echo "Validating and restarting the browser"
    remote "cd $REMOTE_ROOT/website && \
        /opt/morphic/venv/bin/pip install --disable-pip-version-check -q -r requirements-lock.txt && \
        MORPHIC_DATA_PATH=$REMOTE_ROOT/data /opt/morphic/venv/bin/python scripts/validate_data_sources.py && \
        sudo cp deploy/morphic.service /etc/systemd/system/morphic.service && \
        sudo systemctl daemon-reload && \
        sudo systemctl restart morphic && \
        sleep 3 && \
        curl --fail --silent http://127.0.0.1:8501/_stcore/health"

    remote "sudo cp $REMOTE_ROOT/website/deploy/nginx.conf /etc/nginx/sites-available/perturb.dev && \
        sudo nginx -t && sudo systemctl reload nginx"
fi

if [ "$MODE" = "quick" ]; then
    echo
    echo "=== Quick sync complete ==="
    exit 0
fi

echo
echo "[3/4] Syncing complete figure and per-gene analysis collections"
sync_directory "$DATA_ROOT/figures" "$REMOTE_ROOT/data/tf_perturbseq/figures"
mirror_directory \
    "$DATA_ROOT/figures/EBs" \
    "$REMOTE_ROOT/data/tf_perturbseq/figures/EBs" \
    "figures_EBs"
mirror_directory \
    "$DATA_ROOT/results/EBs/differential_expression/deg_tables" \
    "$REMOTE_ROOT/data/tf_perturbseq/results/EBs/differential_expression/deg_tables" \
    "deg_tables"
mirror_directory \
    "$DATA_ROOT/results/EBs/lineage_de" \
    "$REMOTE_ROOT/data/tf_perturbseq/results/EBs/lineage_de" \
    "lineage_de"

echo
echo "[4/4] Syncing stable CellxGene datasets"
sync_cellxgene_file "$EB_CELLXGENE" "$CELLXGENE_DATA_DIR/MORPHIC_EBs_latest.h5ad"
sync_cellxgene_file "$IPSC_CELLXGENE" "$CELLXGENE_DATA_DIR/MORPHIC_iPSC_latest.h5ad"
sync_cellxgene_file "$TIMECOURSE_CELLXGENE" "$CELLXGENE_DATA_DIR/MORPHIC_timecourse_latest.h5ad"

if [ "$DRY_RUN" = "0" ]; then
    remote "find $CELLXGENE_DATA_DIR -maxdepth 1 -name 'MORPHIC_*_latest.h5ad' \
        -printf '%f %s %TY-%Tm-%TdT%TH:%TM:%TSZ\n' | sort; \
        df -h $REMOTE_ROOT $CELLXGENE_DATA_DIR"
fi

echo
echo "=== Full production sync complete ==="
