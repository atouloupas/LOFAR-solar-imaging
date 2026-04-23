from pathlib import Path


def test_expected_directories_exist():
    expected = [
        "configs",
        "data/raw",
        "data/interim",
        "data/processed",
        "docs",
        "outputs/images_fits",
        "outputs/images_png",
        "outputs/movies",
        "outputs/solutions",
        "scripts/legacy",
        "src/lofar_solar_imaging",
        "tests",
    ]
    for path in expected:
        assert Path(path).exists(), f"Missing expected path: {path}"
