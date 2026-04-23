#!/bin/bash
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"

export OPENBLAS_NUM_THREADS="${OPENBLAS_NUM_THREADS:-1}"

usage() {
  cat <<USAGE
Usage: $0 --cal-ms <file> --sun-ms <file> [options]

Options:
  --cal-ms <file>        Calibrator measurement set.
  --sun-ms <file>        Sun (target) measurement set.
  --source-db <file>     Sky model database (default: Cas-sources.txt).
  --solutions-dir <dir>  Output directory for solutions (default: outputs/solutions).
  --skip-predict         Skip predict step.
  --skip-gaincal         Skip gain calibration step.
  --skip-apply           Skip applycal step.
  -h, --help             Show this help.
USAGE
}

CAL_FILE=""
SUN_FILE=""
SOURCE_DB="$REPO_ROOT/Cas-sources.txt"
SOLUTIONS_DIR="$REPO_ROOT/outputs/solutions"
SKIP_PREDICT=0
SKIP_GAINCAL=0
SKIP_APPLY=0

while [[ $# -gt 0 ]]; do
  case "$1" in
    --cal-ms) CAL_FILE="$2"; shift 2 ;;
    --sun-ms) SUN_FILE="$2"; shift 2 ;;
    --source-db) SOURCE_DB="$2"; shift 2 ;;
    --solutions-dir) SOLUTIONS_DIR="$2"; shift 2 ;;
    --skip-predict) SKIP_PREDICT=1; shift ;;
    --skip-gaincal) SKIP_GAINCAL=1; shift ;;
    --skip-apply) SKIP_APPLY=1; shift ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown option: $1"; usage; exit 1 ;;
  esac
done

if [[ -z "$CAL_FILE" || -z "$SUN_FILE" ]]; then
  echo "Error: --cal-ms and --sun-ms are required"
  usage
  exit 1
fi

if [[ ! -f "$SOURCE_DB" ]]; then
  echo "Source database '$SOURCE_DB' not found"
  exit 1
fi

mkdir -p "$SOLUTIONS_DIR"
BASE_NAME="$(basename "$SUN_FILE")"
BASE_NAME="${BASE_NAME%%.*}"
PARMDB_FILE="$SOLUTIONS_DIR/solutions_${BASE_NAME}.h5"

for table in "$CAL_FILE" "$SUN_FILE"; do
  taql "update $table set FLAG_ROW=false"
  taql "update $table set WEIGHT_SPECTRUM=1"
  taql "update $table set FLAG=false"
  taql "update $table set FLAG=true where ANTENNA1=ANTENNA2"
done

if [[ $SKIP_PREDICT -eq 0 ]]; then
  echo "Running predict on $CAL_FILE"
  DP3 steps=[predict] \
    msin="$CAL_FILE" \
    msout=. \
    msout.datacolumn=MODEL_DATA \
    predict.sourcedb="$SOURCE_DB" \
    predict.usebeammodel=true
fi

if [[ $SKIP_GAINCAL -eq 0 ]]; then
  echo "Detecting outlier stations from $SUN_FILE"
  python3 "$REPO_ROOT/src/lofar_solar_imaging/utils/station_outliers.py" "$SUN_FILE"
  mapfile -t stations < "stations/stations_${BASE_NAME}.txt"

  baseline_flags=""
  for station in "${stations[@]}"; do
    baseline_flags+="[$station,*],[*,$station],"
  done
  baseline_flags="${baseline_flags%,}"

  if [[ -n "$baseline_flags" ]]; then
    preflag_baseline="[$baseline_flags,[RS*,RS*]]"
  else
    preflag_baseline="[[RS*,RS*]]"
  fi

  echo "Running gain calibration"
  DP3 steps=[preflagger,gaincal] \
    preflagger.baseline="$preflag_baseline" \
    gaincal.caltype=diagonal \
    gaincal.onebeamperpatch=true \
    gaincal.parmdb="$PARMDB_FILE" \
    gaincal.sourcedb="$SOURCE_DB" \
    gaincal.applysolution=true \
    gaincal.usebeammodel=true \
    gaincal.uvlambdamin=120 \
    msin="$CAL_FILE" \
    msout=. \
    msout.datacolumn=CALIBRATED_DATA
fi

if [[ $SKIP_APPLY -eq 0 ]]; then
  if [[ -f interpolate_solutions.py ]]; then
    python3 interpolate_solutions.py "$PARMDB_FILE" --interpolate --save
  fi

  echo "Applying calibration solutions to $SUN_FILE"
  DP3 steps=[applyphase,applyamp,applybeam] \
    applybeam.updateweights=true \
    applyphase.type=applycal \
    applyphase.interpolation=nearest \
    applyphase.updateweights=true \
    applyphase.parmdb="$PARMDB_FILE" \
    applyphase.correction=phase000 \
    applyphase.soltab=sol000 \
    applyamp.type=applycal \
    applyamp.interpolation=nearest \
    applyamp.updateweights=true \
    applyamp.parmdb="$PARMDB_FILE" \
    applyamp.correction=amplitude000 \
    applyamp.soltab=sol000 \
    msin="$SUN_FILE" \
    msout=. \
    msout.datacolumn=CALIBRATED_DATA
fi

echo "Calibration pipeline completed."
