# LOFAR-solar-imaging

Data analysis and imaging of LOFAR data.

## Project layout

Core Python logic now lives in a package under `src/lofar_solar_imaging/`:

- `imaging/`: imaging-domain utilities such as calibrator model flux helpers.
- `reprojection/`: FITS reprojection to helioprojective coordinates.
- `plotting/`: plotting helpers for helioprojective and FITS overview images.
- `qc/`: quality control helpers (e.g., station outlier detection).
- `utils/`: shared helpers including YAML config loading.

Thin command-line entry points are in `scripts/` and only parse arguments + invoke package functions.

## Configuration

Default runtime configuration lives in `configs/default.yaml`.

The CLIs read this file to resolve:

- runtime storage paths;
- field-of-view defaults;
- QC thresholds and output locations.

Override defaults by passing `--config <path>` to any script in `scripts/`.

## Runtime storage expectations

These directories are reserved for runtime artifacts and are intentionally gitignored:

- `data/raw`: raw ingested data.
- `data/intermediate`: intermediate pipeline outputs.
- `data/products`: final products, reports, and exported artifacts.

The repository tracks only `.gitkeep` placeholders so the directory structure exists.

## Tests

Tests are organized as:

- `tests/unit`: unit tests for package-level logic.
- `tests/integration`: smoke-level integration tests for major pipeline stages.

Run with:

```bash
PYTHONPATH=src pytest
```
