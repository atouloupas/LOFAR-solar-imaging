from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from .common import ensure_exists, run_command, strip_all_suffixes


@dataclass(slots=True)
class ImagingOptions:
    ms: Path
    image_type: str = "calibrated"
    intervals_out: int = 15
    output_dir: Path = Path("outputs/images_fits")
    niter: int = 5000
    repo_root: Path = Path.cwd()


def run_imaging(options: ImagingOptions) -> None:
    ensure_exists(options.ms, description="Input MS")

    data_column_lookup = {
        "raw": "DATA",
        "model": "MODEL_DATA",
        "calibrated": "CALIBRATED_DATA",
    }
    if options.image_type not in data_column_lookup:
        raise ValueError("--type must be one of: raw, model, calibrated")

    base_name = strip_all_suffixes(options.ms)
    file_name = f"{base_name}-{options.image_type}"
    image_dir = options.output_dir / file_name
    image_dir.mkdir(parents=True, exist_ok=True)
    output_name = image_dir / file_name

    run_command(
        [
            "wsclean",
            "-name",
            str(output_name),
            "-mem",
            "50",
            "-weight",
            "briggs",
            "0.2",
            "-scale",
            "30asec",
            "-size",
            "512",
            "512",
            "-pol",
            "I",
            "-multiscale",
            "-data-column",
            data_column_lookup[options.image_type],
            "-niter",
            str(options.niter),
            "-no-reorder",
            "-no-update-model-required",
            "-mgain",
            "0.8",
            "-maxuvw-m",
            "3000",
            "-auto-mask",
            "6",
            "-auto-threshold",
            "2",
            "-fit-beam",
            "-make-psf",
            "-intervals-out",
            str(options.intervals_out),
            str(options.ms),
        ]
    )

    print("WSClean processing completed.")
    if options.intervals_out == 1:
        run_command(["python3", str(options.repo_root / "helioprojective_plot.py"), f"{output_name}-image.fits"])
    else:
        run_command(
            [
                "python3",
                str(options.repo_root / "reproject.py"),
                f"{output_name}-t*-image.fits",
                f"outputs/images_png/{file_name}/",
            ]
        )
