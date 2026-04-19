from lofar_solar_imaging.utils import load_config


def test_load_default_config():
    cfg = load_config('configs/default.yaml')
    assert cfg['paths']['raw'] == 'data/raw'
    assert cfg['imaging']['fov_arcsec'] == 3000
