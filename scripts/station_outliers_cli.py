#!/usr/bin/env python3
"""Thin CLI wrapper for station outlier detection."""

from argparse import ArgumentParser
import numpy as np

from lofar_solar_imaging.qc.station_qc import compute_station_outliers, save_outlier_stations
from lofar_solar_imaging.utils import load_config, ensure_runtime_dirs


def parse_args():
    parser = ArgumentParser(description='Detect station outliers from precomputed flux arrays')
    parser.add_argument('flux_file', help='NumPy .npy file containing station-averaged fluxes')
    parser.add_argument('station_file', help='Text file with station names, one per line')
    parser.add_argument('--config', default='configs/default.yaml')
    parser.add_argument('--std-thresh', default=None, type=float)
    parser.add_argument('--output-file', default=None)
    return parser.parse_args()


def main():
    args = parse_args()
    config = load_config(args.config)
    ensure_runtime_dirs(config)
    with open(args.station_file, 'r') as stream:
        station_names = [line.strip() for line in stream if line.strip()]
    station_fluxes = np.load(args.flux_file)
    std_thresh = args.std_thresh or config['qc']['std_thresh']

    outliers, lower, upper = compute_station_outliers(station_fluxes, std_thresh=std_thresh)
    output_file = args.output_file or config['qc']['outlier_station_file']
    save_outlier_stations(station_names, outliers, output_file)
    print(f'Outliers saved to {output_file} (lower={lower:.6g}, upper={upper:.6g})')


if __name__ == '__main__':
    main()
