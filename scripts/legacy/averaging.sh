#!/bin/bash
set -euo pipefail

export OPENBLAS_NUM_THREADS="${OPENBLAS_NUM_THREADS:-1}"

usage() {
  cat <<USAGE
Usage: $0 --input-dir <dir> [options]

Options:
  --input-dir <dir>      Directory containing *.MS files.
  --sun-sap <SAP>        SAP used for Sun merge output (default: 000).
  --cal-sap <SAP>        SAP used for calibrator merge output (default: 001).
  --obsid <id>           Optional filter for a single observation ID (e.g. L883060).
  --freqres <Hz>         DP3 averager.freqresolution (default: 500000).
  --timeres <sec>        DP3 averager.timeresolution (default: 60).
  --skip-copy            Skip MS -> .CPY stage.
  --skip-average         Skip .CPY -> .AVG stage.
  --skip-merge           Skip SAP merge stage.
  -h, --help             Show this help.
USAGE
}

INPUT_DIR=""
SUN_SAP="000"
CAL_SAP="001"
OBSID=""
FREQRES="500000"
TIMERES="60"
SKIP_COPY=0
SKIP_AVERAGE=0
SKIP_MERGE=0

while [[ $# -gt 0 ]]; do
  case "$1" in
    --input-dir) INPUT_DIR="$2"; shift 2 ;;
    --sun-sap) SUN_SAP="$2"; shift 2 ;;
    --cal-sap) CAL_SAP="$2"; shift 2 ;;
    --obsid) OBSID="$2"; shift 2 ;;
    --freqres) FREQRES="$2"; shift 2 ;;
    --timeres) TIMERES="$2"; shift 2 ;;
    --skip-copy) SKIP_COPY=1; shift ;;
    --skip-average) SKIP_AVERAGE=1; shift ;;
    --skip-merge) SKIP_MERGE=1; shift ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown option: $1"; usage; exit 1 ;;
  esac
done

if [[ -z "$INPUT_DIR" ]]; then
  echo "Error: --input-dir is required"
  usage
  exit 1
fi

shopt -s nullglob
if [[ -n "$OBSID" ]]; then
  MS_FILES=("$INPUT_DIR"/"$OBSID"_SAP*_SB*_uv.MS)
else
  MS_FILES=("$INPUT_DIR"/*.MS)
fi

if [[ ${#MS_FILES[@]} -eq 0 ]]; then
  echo "No measurement sets found in $INPUT_DIR"
  exit 1
fi

if [[ $SKIP_COPY -eq 0 ]]; then
  for ms in "${MS_FILES[@]}"; do
    echo "Copying $ms -> $ms.CPY"
    DP3 msin="$ms" msout="$ms.CPY" steps=[]
  done
fi

if [[ $SKIP_AVERAGE -eq 0 ]]; then
  CPY_FILES=("$INPUT_DIR"/*.MS.CPY)
  for ms in "${CPY_FILES[@]}"; do
    [[ -e "$ms" ]] || continue
    echo "Preparing flags/weights in $ms"
    taql "update $ms set FLAG_ROW=false"
    taql "update $ms set WEIGHT_SPECTRUM=1"
    taql "update $ms set FLAG=false"
    taql "update $ms set FLAG=true where ANTENNA1=ANTENNA2"

    echo "Averaging $ms -> $ms.AVG"
    DP3 msin="$ms" \
      steps=[averager] \
      msout="$ms.AVG" \
      averager.freqresolution="$FREQRES" \
      averager.timeresolution="$TIMERES"
  done
fi

if [[ $SKIP_MERGE -eq 0 ]]; then
  DATE_TAG="$(basename "$INPUT_DIR")"

  SUN_INPUTS=("$INPUT_DIR"/*_SAP"$SUN_SAP"_*.AVG)
  if [[ ${#SUN_INPUTS[@]} -gt 0 ]]; then
    DP3 msin=["$INPUT_DIR"/*_SAP"$SUN_SAP"_*.AVG] steps=[] msout="joined_sun_${DATE_TAG}.MS"
    echo "Created joined_sun_${DATE_TAG}.MS"
  else
    echo "No averaged Sun files found for SAP${SUN_SAP}; skipping Sun merge."
  fi

  CAL_INPUTS=("$INPUT_DIR"/*_SAP"$CAL_SAP"_*.AVG)
  if [[ ${#CAL_INPUTS[@]} -gt 0 ]]; then
    DP3 msin=["$INPUT_DIR"/*_SAP"$CAL_SAP"_*.AVG] steps=[] msout="joined_cal_${DATE_TAG}.MS"
    echo "Created joined_cal_${DATE_TAG}.MS"
  else
    echo "No averaged calibrator files found for SAP${CAL_SAP}; skipping calibrator merge."
  fi
fi

echo "Averaging pipeline completed."
