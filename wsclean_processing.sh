#!/bin/bash
set -euo pipefail

export OPENBLAS_NUM_THREADS="${OPENBLAS_NUM_THREADS:-1}"

usage() {
  cat <<USAGE
Usage: $0 --ms <file> [options]

Options:
  --ms <file>            Input measurement set.
  --type <raw|model|calibrated>  Data column selection (default: calibrated).
  --intervals-out <n>    WSClean intervals-out value (default: 15).
  --output-dir <dir>     Output root directory for FITS images (default: images_fits).
  --niter <n>            Number of clean iterations (default: 5000).
  -h, --help             Show this help.
USAGE
}

FILE=""
TYPE="calibrated"
INTERVALS_OUT=15
OUTPUT_DIR="images_fits"
NITER=5000

while [[ $# -gt 0 ]]; do
  case "$1" in
    --ms) FILE="$2"; shift 2 ;;
    --type) TYPE="$2"; shift 2 ;;
    --intervals-out) INTERVALS_OUT="$2"; shift 2 ;;
    --output-dir) OUTPUT_DIR="$2"; shift 2 ;;
    --niter) NITER="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown option: $1"; usage; exit 1 ;;
  esac
done

if [[ -z "$FILE" ]]; then
  echo "Error: --ms is required"
  usage
  exit 1
fi

BASE_NAME="$(basename "$FILE")"
BASE_NAME="${BASE_NAME%%.*}"
FILE_NAME="${BASE_NAME}-${TYPE}"
mkdir -p "$OUTPUT_DIR/$FILE_NAME"

case "$TYPE" in
  model) DATA_COLUMN="MODEL_DATA" ;;
  calibrated) DATA_COLUMN="CALIBRATED_DATA" ;;
  raw) DATA_COLUMN="DATA" ;;
  *)
    echo "Invalid --type '$TYPE'. Use raw, model, or calibrated."
    exit 1
    ;;
esac

OUTPUT_NAME="$OUTPUT_DIR/$FILE_NAME/$FILE_NAME"

wsclean -name "$OUTPUT_NAME" \
  -mem 50 \
  -weight briggs 0.2 \
  -scale 30asec \
  -size 512 512 \
  -pol I \
  -multiscale \
  -data-column "$DATA_COLUMN" \
  -niter "$NITER" \
  -no-reorder \
  -no-update-model-required \
  -mgain 0.8 \
  -maxuvw-m 3000 \
  -auto-mask 6 \
  -auto-threshold 2 \
  -fit-beam \
  -make-psf \
  -intervals-out "$INTERVALS_OUT" \
  "$FILE"

echo "WSClean processing completed."

if [[ "$INTERVALS_OUT" -eq 1 ]]; then
  python3 helioprojective_plot.py "$OUTPUT_NAME-image.fits"
else
  python3 reproject.py "$OUTPUT_NAME"-t*-image.fits "images_png/$FILE_NAME/"
fi
