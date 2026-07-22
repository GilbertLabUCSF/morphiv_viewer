#!/bin/bash
# Historical alias retained for existing notes and shell history.
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
exec "$SCRIPT_DIR/rsync_to_gcp.sh" --full "$@"
