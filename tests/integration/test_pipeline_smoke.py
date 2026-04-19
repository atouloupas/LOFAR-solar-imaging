import numpy as np

from lofar_solar_imaging.pipeline import imaging_stage, qc_stage


def test_imaging_stage_smoke():
    result = imaging_stage('CasA', np.array([50e6]))
    assert result[0] > 0


def test_qc_stage_smoke():
    outliers, _, _ = qc_stage(np.array([1.0, 1.1, 0.95, 6.0]), std_thresh=1.0)
    assert outliers.any()
