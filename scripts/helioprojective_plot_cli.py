#!/usr/bin/env python3
"""Thin CLI wrapper for helioprojective plotting."""

import os
from argparse import ArgumentParser

from lofar_solar_imaging.plotting import plot_helioprojective_fits
from lofar_solar_imaging.utils import load_config, ensure_runtime_dirs


def parse_args():
    parser = ArgumentParser(description='Generate helioprojective plot from FITS')
    parser.add_argument('filename', help='Input FITS file')
    parser.add_argument('--config', default='configs/default.yaml')
    parser.add_argument('--output', default=None)
    return parser.parse_args()


def main():
    args = parse_args()
    config = load_config(args.config)
    ensure_runtime_dirs(config)
    basename = os.path.splitext(os.path.basename(args.filename))[0]
    output = args.output or os.path.join(config['paths']['products'], f'helioprojective_{basename}.png')
    path = plot_helioprojective_fits(args.filename, output, fov=config['imaging']['fov_arcsec'])
    print(f'Helioprojective image completed: {path}')


if __name__ == '__main__':
    main()
