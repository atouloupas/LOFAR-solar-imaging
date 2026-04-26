from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from .common import ensure_exists, run_command, strip_all_suffixes


@dataclass(slots=True)
class CalibrationOptions:
    cal_ms: Path
    sun_ms: Path
    source_db: Path
    solutions_dir: Path
    skip_predict: bool = False
    skip_gaincal: bool = False
    skip_apply: bool = False


def _prepare_flags(ms_path: Path) -> None:
    run_command(["taql", f"update {ms_path} set FLAG_ROW=false"])
    run_command(["taql", f"update {ms_path} set WEIGHT_SPECTRUM=1"])
    run_command(["taql", f"update {ms_path} set FLAG=false"])
    run_command(["taql", f"update {ms_path} set FLAG=true where ANTENNA1=ANTENNA2"])


def _build_preflag_baseline(base_name: str) -> str:
    stations_file = Path("stations") / f"stations_{base_name}.txt"
    stations: list[str] = []
    if stations_file.exists():
        stations = [line.strip() for line in stations_file.read_text().splitlines() if line.strip()]

    flagged = ",".join(f"[{station},*],[*,{station}]" for station in stations)
    if flagged:
        return f"[{flagged},[RS*,RS*]]"
    return "[[RS*,RS*]]"


def run_calibration(options: CalibrationOptions) -> Path:
    ensure_exists(options.cal_ms, description="Calibrator MS")
    ensure_exists(options.sun_ms, description="Sun MS")
    ensure_exists(options.source_db, description="Source DB")

    options.solutions_dir.mkdir(parents=True, exist_ok=True)
    base_name = strip_all_suffixes(options.sun_ms)
    parmdb_file = options.solutions_dir / f"solutions_{base_name}.h5"

    for table in (options.cal_ms, options.sun_ms):
        _prepare_flags(table)

    if not options.skip_predict:
        run_command(
            [
                "DP3",
                "steps=[predict]",
                f"msin={options.cal_ms}",
                "msout=.",
                "msout.datacolumn=MODEL_DATA",
                f"predict.sourcedb={options.source_db}",
                "predict.usebeammodel=true",
            ]
        )

    if not options.skip_gaincal:
        run_command(["python3", "-m", "lofar_solar_imaging.utils.station_outliers", str(options.sun_ms)])
        preflag_baseline = _build_preflag_baseline(base_name)
        run_command(
            [
                "DP3",
                "steps=[preflagger,gaincal]",
                f"preflagger.baseline={preflag_baseline}",
                "gaincal.caltype=diagonal",
                "gaincal.onebeamperpatch=true",
                f"gaincal.parmdb={parmdb_file}",
                f"gaincal.sourcedb={options.source_db}",
                "gaincal.applysolution=true",
                "gaincal.usebeammodel=true",
                "gaincal.uvlambdamin=120",
                f"msin={options.cal_ms}",
                "msout=.",
                "msout.datacolumn=CALIBRATED_DATA",
            ]
        )

    if not options.skip_apply:
        interpolate_script = Path("interpolate_solutions.py")
        if interpolate_script.exists():
            run_command(["python3", str(interpolate_script), str(parmdb_file), "--interpolate", "--save"])

        run_command(
            [
                "DP3",
                "steps=[applyphase,applyamp,applybeam]",
                "applybeam.updateweights=true",
                "applyphase.type=applycal",
                "applyphase.interpolation=nearest",
                "applyphase.updateweights=true",
                f"applyphase.parmdb={parmdb_file}",
                "applyphase.correction=phase000",
                "applyphase.soltab=sol000",
                "applyamp.type=applycal",
                "applyamp.interpolation=nearest",
                "applyamp.updateweights=true",
                f"applyamp.parmdb={parmdb_file}",
                "applyamp.correction=amplitude000",
                "applyamp.soltab=sol000",
                f"msin={options.sun_ms}",
                "msout=.",
                "msout.datacolumn=CALIBRATED_DATA",
            ]
        )

    print("Calibration pipeline completed.")
    return parmdb_file
