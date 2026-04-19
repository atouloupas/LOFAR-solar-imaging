"""Small orchestration helpers for pipeline stages."""

from lofar_solar_imaging.imaging import model_flux
from lofar_solar_imaging.qc import compute_station_outliers


def imaging_stage(calibrator: str, frequency_hz):
    return model_flux(calibrator, frequency_hz)


def qc_stage(station_fluxes, std_thresh: float):
    return compute_station_outliers(station_fluxes, std_thresh=std_thresh)
