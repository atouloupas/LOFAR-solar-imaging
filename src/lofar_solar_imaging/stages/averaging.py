from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from .common import as_dp3_msin, run_command


@dataclass(slots=True)
class AveragingOptions:
    input_dir: Path
    sun_sap: str = "000"
    cal_sap: str = "001"
    obsid: str | None = None
    freqres: str = "500000"
    timeres: str = "60"
    skip_copy: bool = False
    skip_average: bool = False
    skip_merge: bool = False


def _find_ms_files(options: AveragingOptions) -> list[Path]:
    pattern = f"{options.obsid}_SAP*_SB*_uv.MS" if options.obsid else "*.MS"
    files = sorted(options.input_dir.glob(pattern))
    if not files:
        raise FileNotFoundError(f"No measurement sets found in {options.input_dir}")
    return files


def _prepare_flags(ms_path: Path) -> None:
    run_command(["taql", f"update {ms_path} set FLAG_ROW=false"])
    run_command(["taql", f"update {ms_path} set WEIGHT_SPECTRUM=1"])
    run_command(["taql", f"update {ms_path} set FLAG=false"])
    run_command(["taql", f"update {ms_path} set FLAG=true where ANTENNA1=ANTENNA2"])


def run_averaging(options: AveragingOptions) -> None:
    ms_files = _find_ms_files(options)

    if not options.skip_copy:
        for ms in ms_files:
            run_command(["DP3", f"msin={ms}", f"msout={ms}.CPY", "steps=[]"])

    if not options.skip_average:
        for copied_ms in sorted(options.input_dir.glob("*.MS.CPY")):
            _prepare_flags(copied_ms)
            run_command(
                [
                    "DP3",
                    f"msin={copied_ms}",
                    "steps=[averager]",
                    f"msout={copied_ms}.AVG",
                    f"averager.freqresolution={options.freqres}",
                    f"averager.timeresolution={options.timeres}",
                ]
            )

    if not options.skip_merge:
        date_tag = options.input_dir.name

        sun_inputs = sorted(options.input_dir.glob(f"*_SAP{options.sun_sap}_*.AVG"))
        if sun_inputs:
            run_command(
                [
                    "DP3",
                    f"msin={as_dp3_msin(sun_inputs)}",
                    "steps=[]",
                    f"msout=joined_sun_{date_tag}.MS",
                ]
            )
            print(f"Created joined_sun_{date_tag}.MS")
        else:
            print(f"No averaged Sun files found for SAP{options.sun_sap}; skipping Sun merge.")

        cal_inputs = sorted(options.input_dir.glob(f"*_SAP{options.cal_sap}_*.AVG"))
        if cal_inputs:
            run_command(
                [
                    "DP3",
                    f"msin={as_dp3_msin(cal_inputs)}",
                    "steps=[]",
                    f"msout=joined_cal_{date_tag}.MS",
                ]
            )
            print(f"Created joined_cal_{date_tag}.MS")
        else:
            print(f"No averaged calibrator files found for SAP{options.cal_sap}; skipping calibrator merge.")

    print("Averaging pipeline completed.")
