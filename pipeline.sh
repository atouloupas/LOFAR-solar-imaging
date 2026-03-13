#!/bin/bash
set -euo pipefail

usage() {
  cat <<USAGE
Usage: $0 [steps] --input-dir <dir> [selection options]

Step flags (at least one required):
  --average              Run averaging.sh.
  --calibrate            Run calibration.sh.
  --image                Run wsclean_processing.sh.

Selection/input options:
  --input-dir <dir>      Directory containing MS files.
  --sun-sap <SAP>        SAP ID for Sun/target (default: 000).
  --cal-sap <SAP>        SAP ID for calibrator (default: 001).
  --sb <SB>              Optional SB filter (e.g. 230) for single-file calibration/imaging.
  --obsid <Lxxxxxx>      Optional OBSID filter.
  --sun-ms <file>        Override auto-selected Sun MS for calibration/imaging.
  --cal-ms <file>        Override auto-selected calibrator MS for calibration.

Pass-through options:
  --source-db <file>     Source DB path for calibration (default: Cas-sources.txt).
  --solutions-dir <dir>  Solution directory (default: ./solutions).
  --image-type <type>    raw|model|calibrated for imaging (default: calibrated).
  --intervals-out <n>    intervals-out for imaging (default: 15).

Examples:
  $0 --average --input-dir HBA_files/20230311
  $0 --calibrate --image --input-dir HBA_files/20230311 --sb 230 --sun-sap 000 --cal-sap 001
USAGE
}

RUN_AVERAGE=0
RUN_CALIBRATE=0
RUN_IMAGE=0
INPUT_DIR=""
SUN_SAP="000"
CAL_SAP="001"
SB=""
OBSID=""
SUN_MS=""
CAL_MS=""
SOURCE_DB="Cas-sources.txt"
SOLUTIONS_DIR="./solutions"
IMAGE_TYPE="calibrated"
INTERVALS_OUT="15"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --average) RUN_AVERAGE=1; shift ;;
    --calibrate) RUN_CALIBRATE=1; shift ;;
    --image) RUN_IMAGE=1; shift ;;
    --input-dir) INPUT_DIR="$2"; shift 2 ;;
    --sun-sap) SUN_SAP="$2"; shift 2 ;;
    --cal-sap) CAL_SAP="$2"; shift 2 ;;
    --sb) SB="$2"; shift 2 ;;
    --obsid) OBSID="$2"; shift 2 ;;
    --sun-ms) SUN_MS="$2"; shift 2 ;;
    --cal-ms) CAL_MS="$2"; shift 2 ;;
    --source-db) SOURCE_DB="$2"; shift 2 ;;
    --solutions-dir) SOLUTIONS_DIR="$2"; shift 2 ;;
    --image-type) IMAGE_TYPE="$2"; shift 2 ;;
    --intervals-out) INTERVALS_OUT="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown option: $1"; usage; exit 1 ;;
  esac
done

if [[ $RUN_AVERAGE -eq 0 && $RUN_CALIBRATE -eq 0 && $RUN_IMAGE -eq 0 ]]; then
  echo "Error: select at least one step (--average, --calibrate, --image)."
  usage
  exit 1
fi

if [[ -z "$INPUT_DIR" ]]; then
  echo "Error: --input-dir is required"
  exit 1
fi

resolve_ms() {
  local sap="$1"
  local sb_pattern="*"
  local obs_pattern="*"

  if [[ -n "$SB" ]]; then
    sb_pattern="SB${SB}_uv.MS.CPY.AVG"
  else
    sb_pattern="SB*_uv.MS.CPY.AVG"
  fi

  if [[ -n "$OBSID" ]]; then
    obs_pattern="$OBSID"
  fi

  local matches=("$INPUT_DIR"/"$obs_pattern"_SAP"$sap"_"$sb_pattern")
  if [[ ${#matches[@]} -eq 0 || ! -e "${matches[0]}" ]]; then
    return 1
  fi
  printf '%s\n' "${matches[0]}"
}

if [[ $RUN_AVERAGE -eq 1 ]]; then
  avg_cmd=(./averaging.sh --input-dir "$INPUT_DIR" --sun-sap "$SUN_SAP" --cal-sap "$CAL_SAP")
  if [[ -n "$OBSID" ]]; then
    avg_cmd+=(--obsid "$OBSID")
  fi
  "${avg_cmd[@]}"
fi

if [[ -z "$SUN_MS" ]]; then
  if SUN_CANDIDATE="$(resolve_ms "$SUN_SAP")"; then
    SUN_MS="$SUN_CANDIDATE"
  fi
fi

if [[ -z "$CAL_MS" ]]; then
  if CAL_CANDIDATE="$(resolve_ms "$CAL_SAP")"; then
    CAL_MS="$CAL_CANDIDATE"
  fi
fi

if [[ $RUN_CALIBRATE -eq 1 ]]; then
  if [[ -z "$SUN_MS" || -z "$CAL_MS" ]]; then
    echo "Error: could not auto-resolve SUN/CAL MS. Pass --sun-ms and --cal-ms explicitly."
    exit 1
  fi
  ./calibration.sh --cal-ms "$CAL_MS" --sun-ms "$SUN_MS" --source-db "$SOURCE_DB" --solutions-dir "$SOLUTIONS_DIR"
fi

if [[ $RUN_IMAGE -eq 1 ]]; then
  if [[ -z "$SUN_MS" ]]; then
    echo "Error: could not auto-resolve Sun MS for imaging. Pass --sun-ms explicitly."
    exit 1
  fi
  ./wsclean_processing.sh --ms "$SUN_MS" --type "$IMAGE_TYPE" --intervals-out "$INTERVALS_OUT"
fi

echo "Requested pipeline steps completed."
