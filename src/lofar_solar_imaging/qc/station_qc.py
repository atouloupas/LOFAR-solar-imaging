from pathlib import Path
from typing import Iterable

import numpy as np
from astropy.time import Time


def mjd_to_utc(time_values, str_format):
    """Convert MJD-like seconds values to UTC strings."""
    time_mjd = np.asarray(time_values) / (24 * 3600)
    time_array = Time(time_mjd, format='mjd')
    return time_array.strftime(str_format)


def compute_station_outliers(station_fluxes: np.ndarray, std_thresh: float = 2.0):
    """Return outlier mask and thresholds based on standard deviation."""
    avg = float(np.mean(station_fluxes))
    std = float(np.std(station_fluxes))
    upper = avg + std_thresh * std
    lower = avg - std_thresh * std
    outliers = (station_fluxes > upper) | (station_fluxes < lower)
    return outliers, lower, upper


def save_outlier_stations(station_names: Iterable[str], outlier_mask, output_file: str):
    """Persist outlier station names to text file."""
    output = Path(output_file)
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open('w') as stream:
        for name, is_outlier in zip(station_names, outlier_mask):
            if is_outlier:
                stream.write(f"{name}\n")
