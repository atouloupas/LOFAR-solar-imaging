import numpy as np

from lofar_solar_imaging.qc import compute_station_outliers


def test_compute_station_outliers_detects_extreme_values():
    data = np.array([1.0, 1.1, 0.9, 10.0])
    outliers, lower, upper = compute_station_outliers(data, std_thresh=1.0)
    assert outliers.tolist() == [False, False, False, True]
    assert lower < upper
