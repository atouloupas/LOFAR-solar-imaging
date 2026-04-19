#!/usr/bin/env python3
"""Thin CLI wrapper for reprojection pipeline."""

from argparse import ArgumentParser
import os
from astropy.io import fits

from lofar_solar_imaging.reprojection import reproject_to_heliocentric_frame
from lofar_solar_imaging.utils import load_config, ensure_runtime_dirs


def parse_args():
    parser = ArgumentParser(description='Reproject FITS image to helioprojective coordinates')
    parser.add_argument('input', nargs='+', help='Input FITS files')
    parser.add_argument('--config', default='configs/default.yaml')
    parser.add_argument('--output-dir', default=None)
    return parser.parse_args()


def main():
    args = parse_args()
    config = load_config(args.config)
    ensure_runtime_dirs(config)
    output_dir = args.output_dir or config['paths']['intermediate']
    os.makedirs(output_dir, exist_ok=True)

    for input_path in sorted(args.input):
        with fits.open(input_path) as hdul:
            header, data = hdul[0].header, hdul[0].data
        solar_map = reproject_to_heliocentric_frame(header, data, fov=config['imaging']['fov_arcsec'])
        output = os.path.join(output_dir, os.path.basename(input_path).replace('.fits', '-reprojected.fits'))
        solar_map.save(output, overwrite=True)
        print(f'Wrote {output}')


if __name__ == '__main__':
    main()
