from __future__ import annotations

import argparse
from pathlib import Path

from .averaging import AveragingOptions, run_averaging
from .calibration import CalibrationOptions, run_calibration
from .imaging import ImagingOptions, run_imaging


def _resolve_ms(input_dir: Path, sap: str, sb: str | None, obsid: str | None) -> Path | None:
    sb_pattern = f"SB{sb}_uv.MS.CPY.AVG" if sb else "SB*_uv.MS.CPY.AVG"
    obs_pattern = obsid if obsid else "*"
    matches = sorted(input_dir.glob(f"{obs_pattern}_SAP{sap}_{sb_pattern}"))
    return matches[0] if matches else None


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="LOFAR solar imaging pipeline")

    parser.add_argument("--average", action="store_true", help="Run averaging stage")
    parser.add_argument("--calibrate", action="store_true", help="Run calibration stage")
    parser.add_argument("--image", action="store_true", help="Run imaging stage")

    parser.add_argument("--input-dir", required=True, type=Path, help="Directory containing MS files")
    parser.add_argument("--sun-sap", default="000")
    parser.add_argument("--cal-sap", default="001")
    parser.add_argument("--sb")
    parser.add_argument("--obsid")
    parser.add_argument("--sun-ms", type=Path)
    parser.add_argument("--cal-ms", type=Path)

    parser.add_argument("--source-db", type=Path, default=Path("Cas-sources.txt"))
    parser.add_argument("--solutions-dir", type=Path, default=Path("outputs/solutions"))

    parser.add_argument("--image-type", default="calibrated")
    parser.add_argument("--intervals-out", type=int, default=15)
    parser.add_argument("--output-dir", type=Path, default=Path("outputs/images_fits"))
    parser.add_argument("--niter", type=int, default=5000)

    parser.add_argument("--freqres", default="500000")
    parser.add_argument("--timeres", default="60")
    parser.add_argument("--skip-copy", action="store_true")
    parser.add_argument("--skip-average", action="store_true")
    parser.add_argument("--skip-merge", action="store_true")
    parser.add_argument("--skip-predict", action="store_true")
    parser.add_argument("--skip-gaincal", action="store_true")
    parser.add_argument("--skip-apply", action="store_true")

    return parser


def main() -> None:
    parser = build_parser()
    args = parser.parse_args()

    if not any([args.average, args.calibrate, args.image]):
        parser.error("Select at least one stage: --average, --calibrate, or --image")

    input_dir = args.input_dir
    sun_ms = args.sun_ms or _resolve_ms(input_dir, args.sun_sap, args.sb, args.obsid)
    cal_ms = args.cal_ms or _resolve_ms(input_dir, args.cal_sap, args.sb, args.obsid)

    if args.average:
        run_averaging(
            AveragingOptions(
                input_dir=input_dir,
                sun_sap=args.sun_sap,
                cal_sap=args.cal_sap,
                obsid=args.obsid,
                freqres=args.freqres,
                timeres=args.timeres,
                skip_copy=args.skip_copy,
                skip_average=args.skip_average,
                skip_merge=args.skip_merge,
            )
        )
        sun_ms = sun_ms or _resolve_ms(input_dir, args.sun_sap, args.sb, args.obsid)
        cal_ms = cal_ms or _resolve_ms(input_dir, args.cal_sap, args.sb, args.obsid)

    if args.calibrate:
        if not sun_ms or not cal_ms:
            parser.error("Could not resolve SUN/CAL MS. Pass --sun-ms and --cal-ms explicitly.")
        run_calibration(
            CalibrationOptions(
                cal_ms=cal_ms,
                sun_ms=sun_ms,
                source_db=args.source_db,
                solutions_dir=args.solutions_dir,
                skip_predict=args.skip_predict,
                skip_gaincal=args.skip_gaincal,
                skip_apply=args.skip_apply,
            )
        )

    if args.image:
        if not sun_ms:
            parser.error("Could not resolve SUN MS for imaging. Pass --sun-ms explicitly.")
        run_imaging(
            ImagingOptions(
                ms=sun_ms,
                image_type=args.image_type,
                intervals_out=args.intervals_out,
                output_dir=args.output_dir,
                niter=args.niter,
                repo_root=Path(__file__).resolve().parents[3],
            )
        )

    print("Requested pipeline stages completed.")


if __name__ == "__main__":
    main()
