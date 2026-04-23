# Pipeline architecture

The pipeline is currently orchestrated with shell scripts in `scripts/legacy` and
is progressively being migrated into Python modules under `src/lofar_solar_imaging`.

## Current flow

1. Averaging (`averaging.sh`)
2. Calibration (`calibration.sh`)
3. Imaging (`wsclean_processing.sh`)

`pipeline.sh`/`scripts/run_pipeline.sh` orchestrates these steps.
