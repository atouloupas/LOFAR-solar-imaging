from pathlib import Path
import yaml


def load_config(config_path: str = 'configs/default.yaml'):
    with open(config_path, 'r') as stream:
        return yaml.safe_load(stream)


def ensure_runtime_dirs(config: dict):
    for folder in config['paths'].values():
        Path(folder).mkdir(parents=True, exist_ok=True)
