#!/usr/bin/env python3
"""Thin CLI wrapper for FITS image plotting."""

import os
from argparse import ArgumentParser

from lofar_solar_imaging.plotting import plot_fits_image
from lofar_solar_imaging.utils import load_config, ensure_runtime_dirs


def parse_args():
    parser = ArgumentParser(description='Plot FITS images')
    parser.add_argument('filename', help='FITS file to plot')
    parser.add_argument('--config', default='configs/default.yaml')
    parser.add_argument('--no-cbar', action='store_true')
    parser.add_argument('--output', default=None)
    return parser.parse_args()


def main():
    args = parse_args()
    config = load_config(args.config)
    ensure_runtime_dirs(config)
    basename = os.path.splitext(os.path.basename(args.filename))[0]
    output = args.output or os.path.join(config['paths']['products'], f'{basename}.png')
    plot_fits_image(args.filename, output, with_colorbar=not args.no_cbar)
    print(f'Plot completed: {output}')


if __name__ == '__main__':
    main()
