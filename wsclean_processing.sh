#!/bin/bash
set -euo pipefail

exec "$(dirname "$0")/scripts/legacy/wsclean_processing.sh" "$@"
