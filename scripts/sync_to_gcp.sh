#!/bin/bash
# Backward-compatible entry point for a quick code/core-data deployment.
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
exec "$SCRIPT_DIR/rsync_to_gcp.sh" --quick "$@"
