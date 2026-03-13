# LOFAR-solar-imaging

Data analysis and imaging of LOFAR data.

## Unified pipeline

Use `pipeline.sh` to run averaging, calibration, and imaging in one command while keeping each step optional.

### Step selection

- `--average`: copy MS files, average them, and merge by SAP.
- `--calibrate`: predict sky model, solve gains, and apply solutions.
- `--image`: run WSClean imaging.

You can run any subset of steps by combining flags.

### Example

```bash
./pipeline.sh \
  --average --calibrate --image \
  --input-dir HBA_files/20230311 \
  --sun-sap 000 \
  --cal-sap 001 \
  --sb 230
```

### Individual scripts

All scripts now support CLI arguments:

- `averaging.sh --help`
- `calibration.sh --help`
- `wsclean_processing.sh --help`

## Input naming convention and parsing

Expected MS filename format:

`Laaaaaa_SAPbbb_SBccc_uv.MS`

Example: `L883060_SAP000_SB230_uv.MS`

### Suggested improvements

1. Keep all files for one run in a single date folder (already done).
2. Always pass `--sun-sap`, `--cal-sap`, and `--sb` to avoid accidental file selection.
3. For reproducibility, consider adding a small manifest file (CSV/YAML) with columns:
   - `obsid`
   - `sap`
   - `sb`
   - `role` (`sun` or `calibrator`)
   - `path`

A manifest removes ambiguity when multiple observation IDs/SBs exist in the same directory.
