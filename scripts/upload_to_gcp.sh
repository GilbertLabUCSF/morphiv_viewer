#!/bin/bash
# Backward-compatible entry point for a complete production refresh.
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
exec "$SCRIPT_DIR/rsync_to_gcp.sh" --full "$@"
