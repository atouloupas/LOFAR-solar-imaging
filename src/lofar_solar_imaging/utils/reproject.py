#!/usr/bin/env python3
from astropy.io import fits
from argparse import ArgumentParser
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord, EarthLocation
import sunpy.map
import sunpy.coordinates
import matplotlib.colors as mpl_colors
import matplotlib.patches as mpl_patches
import numpy as np
import astropy.units as u
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.animation as animation
import os.path
from astropy.time import Time
import re
import json

LOFAR_CORE = EarthLocation.from_geocentric(x=3826560.737625 * u.m, y=461033.24679167 * u.m, z=5064904.97320833 * u.m)


def parse_args():
    parser = ArgumentParser(description='Reproject to sun coordinate')
    parser.add_argument('input', help='input file', nargs='+')
    parser.add_argument('output', help='output file')

    return parser.parse_args()


def reproject_to_heliocentric_frame(header, data, fov=3000):
    """
    Reproject the image generated using equatorial coordinates (ICRS) into
    heliocentric coordinates and changes the fix header accordingly.

    Inspired from Laura Hayes gist
    https://gist.github.com/hayesla/a842afd90aa5a01b8b4d59be8b4ceb8c
    """

    obstime = header['date-obs']
    freq = header['crval3'] * u.Hz
    cdelt1, cdelt2 = abs(header['cdelt1']), abs(header['cdelt2'])
    sun_distance = sunpy.coordinates.sun.earth_distance(obstime)
    reference_pixel = u.Quantity([header['crpix1'], header['crpix2']] * u.pixel)
    data = data.squeeze()
    scale = u.Quantity([cdelt1, cdelt2] * u.deg / u.pix)
    lofar_coord = SkyCoord(LOFAR_CORE.get_itrs(Time(obstime)))
    reference_coordinate = SkyCoord(header['crval1'] * u.deg, header['crval2'] * u.deg,
                                    frame='gcrs',
                                    observer=lofar_coord,
                                    obstime=obstime,
                                    distance=sun_distance)
    reference_coordinate_helio = reference_coordinate.transform_to(sunpy.coordinates.frames.Helioprojective)
    P1 = sunpy.coordinates.sun.P(obstime)
    header_heliocentric_frame = sunpy.map.make_fitswcs_header(data, reference_coordinate_helio,
                                                              reference_pixel=reference_pixel,
                                                              scale=scale,
                                                              rotation_angle=-P1,
                                                              wavelength=freq.to(u.MHz).round(2),
                                                              observatory='LOFAR')

    header_heliocentric_frame['HISTORY'] = header['HISTORY']
    header_heliocentric_frame['BMAJ'] = header['BMAJ']
    header_heliocentric_frame['BMIN'] = header['BMIN']
    header_heliocentric_frame['BPA'] = header['BPA']
    b_maj = header['BMAJ']
    b_min = header['BMIN']
    area = (b_maj) * (b_min) * np.pi
    data *= area 
    data /= 1.e4
    image_map = sunpy.map.Map(data, header_heliocentric_frame)
    image_map = image_map.rotate()

    bl = SkyCoord(-fov * u.arcsec, -fov * u.arcsec, frame=image_map.coordinate_frame)
    tr = SkyCoord(fov * u.arcsec, fov * u.arcsec, frame=image_map.coordinate_frame)
    image_map = image_map.submap(bottom_left=bl, top_right=tr)

    return image_map


def image_plotter(fig: plt.Figure, ax: plt.Axes, sunpy_map: sunpy.map.Map):
    fig.set_dpi(200)
    p, *_ = sunpy_map.draw_limb()
    solar_angle = sunpy.coordinates.sun.P(sunpy_map.date)
    scale = 3600 / 6000
    beam0 = mpl_patches.Ellipse((.7, -.7), sunpy_map.meta['bmaj'] * scale, sunpy_map.meta['bmin'] * scale,
                            angle=-(sunpy_map.meta['bpa']*u.deg + 90*u.deg).value, color='w', transform=ax.get_transform('world'))
    ax.add_patch(beam0)
    return [p]


def plot_solar_map(image_map, path_out):
    mapsmax = 0
    mapsmin = 0
    percentiles_max = []
    percentiles_min = []
    for map in image_map:
        percentile = np.percentile(map.data, q=99)
        percentiles_max.append(percentile)
        mapsmax = max([mapsmax, map.max()])

        percentile = np.percentile(map.data, q=1)
        percentiles_min.append(percentile)
        mapsmin = min([mapsmin, map.min()])
    mapsmax = np.percentile(percentiles_max, q=70)
    # mapsmin = np.percentile(percentiles_min, q=5)
    mapsmin = 1.e-3

    for map in image_map:
        map.plot_settings['cmap'] = plt.get_cmap('inferno')
        map.plot_settings['norm'] = mpl_colors.Normalize(vmin=mapsmin, vmax=mapsmax)
        # map.plot_settings['norm'] = mpl_colors.Normalize(vmin=map.min(), vmax=map.max())
    
    ani = image_map.plot(plot_function=image_plotter)
    plt.colorbar(label='Flux [SFU/pixel]')
    Writer = animation.writers['ffmpeg']
    writer = Writer(fps=2, bitrate=1800)
    ani.save(path_out, writer=writer)

def load_fits(path):
    hdu, *_ = fits.open(path)
    return hdu.header, hdu.data


def get_wcs_from_hdu(hdu):
    return WCS(hdu)

FNAME_PATTERN = r'L(?P<sas_id>\d+)_SAP(?P<sap>\d+)_SB(?P<subband>\d+)_uv.*?-t(?P<sample>\d+)'

def extract_metadata(header, fname):
    metadata = {}
    sas_id, sap, subband, sample = re.search(FNAME_PATTERN, fname).groups()
    metadata['sas_id'] = int(sas_id)
    metadata['sap'] = int(sap)
    metadata['subband'] = int(subband)
    metadata['time_sample'] = int(sample)
    metadata["datetime"] = header['DATE-OBS']
    metadata["date"], metadata['time'] = header['DATE-OBS'].split('T')
    metadata["telescope"] = header["TELESCOP"].strip()
    metadata["bmin"] = header["BMIN"]
    metadata["bmax"] = header["BMAJ"]
    metadata["bpa"] = header["BPA"]
    metadata["frequency"] = header["CRVAL3"]
    metadata["target"] = header["OBJECT"]
    return metadata
    
def metadata_for_video(metadatas):
    metadata = dict(metadatas[0])
    metadata['ntimes'] = len(metadatas)
    metadata['start_time'] = metadata.pop('time')
    metadata['start_date'] = metadata.pop('date')
    metadata['start_datetime'] = metadata.pop('datetime')
    metadata.pop('time_sample')
    metadata['end_time'] = metadatas[-1]['time']
    metadata['end_date'] = metadatas[-1]['date']
    metadata['end_datetime'] = metadatas[-1]['datetime']
    metadata['type'] = 'movie'
    return metadata
    

def main():
    args = parse_args()
    solar_maps = []
    solar_meta = []
    full_path = os.path.join(os.getcwd(), args.output)
    os.makedirs(full_path, exist_ok=True)
    inputs = sorted(args.input)
    file_basename = os.path.basename(inputs[0]).split('.')[0].replace("t0000-image", "")
    pattern = file_basename + 't{index:04}-reprojected.fits'
    for index, input_path in enumerate(inputs):
        print('Processing image ', index + 1, 'of', len(args.input))
        output_basename = pattern.format(index=index)
        input_basename = os.path.basename(input_path)
        header, data = load_fits(input_path)
        solar_map = reproject_to_heliocentric_frame(header, data)
        solar_maps += [solar_map]
        output_solar_map_metadata = extract_metadata(header, output_basename)
        input_solar_map_metadata = extract_metadata(header, input_basename)
        solar_meta += [output_solar_map_metadata]
        with open(os.path.join(full_path, output_basename).replace('.fits', '.json'), 'w') as f_stream:
            output_solar_map_metadata['type'] = 'reprojected-image'
            json.dump(output_solar_map_metadata, f_stream)
        with open(os.path.join(full_path, input_basename).replace('.fits', '.json'), 'w') as f_stream:
            output_solar_map_metadata['type'] = 'image'
            json.dump(input_solar_map_metadata, f_stream)
    solar_map_sequence = sunpy.map.Map(solar_maps, sequence=True)
    
    
    solar_map_sequence.save(os.path.join(args.output, pattern), overwrite=True)
    movie_path = os.path.join(args.output, "MOVIE" + file_basename.rstrip('-') + '.mp4')
    with open(movie_path.replace('.mp4', '.json'), 'w') as f_stream:
        movie_metadata = metadata_for_video(solar_meta)
        json.dump(movie_metadata, f_stream)
    plot_solar_map(solar_map_sequence, movie_path)
    print(f"Movie completed: {movie_path}")


if __name__ == '__main__':
    main()
