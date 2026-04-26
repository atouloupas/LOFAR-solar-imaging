# LOFAR-solar-imaging

Data analysis and imaging of LOFAR solar data.

## Repository structure

```text
LOFAR-solar-imaging/
├── configs/                    # Pipeline configuration files
├── data/                       # Suggested data layout (raw/interim/processed)
├── docs/                       # Architecture/runbook docs
├── outputs/                    # Output products (FITS/PNG/movies/solutions)
├── scripts/
│   ├── run_pipeline.sh         # Preferred shell entrypoint
│   └── legacy/                 # Stage shell scripts
├── src/lofar_solar_imaging/
│   ├── plotting/               # Plotting utilities
│   ├── stages/                 # Stage modules (reserved for migration)
│   └── utils/                  # Processing/calibration helpers
└── tests/
```

## Pipeline entrypoints

The pipeline now runs through a Python entrypoint (`lofar_solar_imaging.stages.pipeline`) and keeps shell wrappers for convenience:

- `./pipeline.sh`
- `./scripts/run_pipeline.sh`
- `python3 -m lofar_solar_imaging.stages.pipeline`

Legacy shell implementations remain under `scripts/legacy/` for reference and fallback.

## Unified pipeline

Use `pipeline.sh` to run averaging, calibration, and imaging in one command while keeping each step optional.

### Step selection

- `--average`: copy MS files, average them, and merge by SAP.
- `--calibrate`: predict sky model, solve gains, and apply solutions.
- `--image`: run WSClean imaging.

You can run any subset of steps by combining flags.

### Example

```bash
./scripts/run_pipeline.sh \
  --average --calibrate --image \
  --input-dir data/raw/20230311 \
  --sun-sap 000 \
  --cal-sap 001 \
  --sb 230
```

## Legacy compatibility wrappers

Top-level wrappers are kept for backwards compatibility:

- `averaging.sh`
- `calibration.sh`
- `wsclean_processing.sh`
- `fits_plot_all.py`
- `helioprojective_plot.py`
- `reproject.py`
- `inspect_ms.py`
- `station_outliers.py`

## Input naming convention and parsing

Expected MS filename format:

`Laaaaaa_SAPbbb_SBccc_uv.MS`

Example: `L883060_SAP000_SB230_uv.MS`
