#!/usr/bin/env python3
import runpy
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent / "src"))
runpy.run_module("lofar_solar_imaging.utils.reproject", run_name="__main__")
