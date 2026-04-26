# Pipeline architecture

The pipeline is orchestrated by Python modules under `src/lofar_solar_imaging/stages`
with shell wrappers kept as lightweight entrypoints.

## Current flow

1. Averaging (`lofar_solar_imaging.stages.averaging`)
2. Calibration (`lofar_solar_imaging.stages.calibration`)
3. Imaging (`lofar_solar_imaging.stages.imaging`)

`pipeline.sh`/`scripts/run_pipeline.sh` call `python3 -m lofar_solar_imaging.stages.pipeline`.
Legacy shell scripts remain in `scripts/legacy/` as compatibility fallback.
