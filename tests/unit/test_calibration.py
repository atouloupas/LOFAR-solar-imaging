import numpy as np

from lofar_solar_imaging.imaging import model_flux


def test_model_flux_returns_positive_values():
    values = model_flux('CasA', np.array([30e6, 60e6]))
    assert values.shape == (2,)
    assert np.all(values > 0)
