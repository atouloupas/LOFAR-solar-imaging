import sunpy.map
from sunpy.coordinates import frames, sun
from astropy.coordinates import SkyCoord, EarthLocation
from astropy import units as u
from astropy.io import fits
from astropy.time import Time
import matplotlib.pyplot as plt
import matplotlib.patches as mpl_patches
import matplotlib.colors as mpl_colors
import numpy as np


def model_flux(calibrator, frequency):
    '''
    Calculates the model matrix for flux calibration for a range of known calibrators:
    J0133-3629, 3C48, Fornax A, 3C 123, J0444+2809, 3C138, Pictor A, Taurus A, 3C147, 3C196, Hydra A, Virgo A,
    3C286, 3C295, Hercules A, 3C353, 3C380, Cygnus A, 3C444, Cassiopeia A

    Input: the calibrator name, frequency range, and time range
    Output: the calibration matrix (in sfu)

    source https://arxiv.org/pdf/1609.05940.pdf

    credits:
     Peijin Zhang, Pietro Zucca
    '''
    parameters = []

    cal_dict = {'J0133-3629': [1.0440, -0.662, -0.225],
                '3C48': [1.3253, -0.7553, -0.1914, 0.0498],
                'ForA': [2.218, -0.661],
                '3C123': [1.8017, -0.7884, -0.1035, -0.0248, 0.0090],
                'J0444-2809': [0.9710, -0.894, -0.118],
                '3C138': [1.0088, -0.4981, -0.155, -0.010, 0.022, ],
                'PicA': [1.9380, -0.7470, -0.074],
                'TauA': [2.9516, -0.217, -0.047, -0.067],
                '3C247': [1.4516, -0.6961, -0.201, 0.064, -0.046, 0.029],
                '3C196': [1.2872, -0.8530, -0.153, -0.0200, 0.0201],
                'HydA': [1.7795, -0.9176, -0.084, -0.0139, 0.030],
                'VirA': [2.4466, -0.8116, -0.048],
                '3C286': [1.2481, -0.4507, -0.1798, 0.0357],
                '3C295': [1.4701, -0.7658, -0.2780, -0.0347, 0.0399],
                'HerA': [1.8298, -1.0247, -0.0951],
                '3C353': [1.8627, -0.6938, -0.100, -0.032],
                '3C380': [1.2320, -0.791, 0.095, 0.098, -0.18, -0.16],
                '3C444': [3.3498, -1.0022, -0.22, 0.023, 0.043],
                'CasA': [3.3584, -0.7518, -0.035, -0.071]}
    if calibrator in cal_dict.keys():
        parameters = cal_dict[calibrator]
    else:
        raise ValueError(calibrator, "is not in the calibrators list")

    flux_model = 0
    frequency = frequency / 10 ** 9  # convert from Hz to GHz
    for j, p in enumerate(parameters):
        flux_model += p * np.log10(frequency) ** j
    flux_model = 10 ** flux_model  # because at first the flux is in log10
    return flux_model * 10 ** (-4)  # convert form Jy to sfu


LOFAR_CORE = EarthLocation.from_geocentric(x=3826560.737625 * u.m, y=461033.24679167 * u.m, z=5064904.97320833 * u.m)

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
    sun_distance = sun.earth_distance(obstime)
    reference_pixel = u.Quantity([header['crpix1'], header['crpix2']] * u.pixel)
    data = data.squeeze()
    scale = u.Quantity([cdelt1, cdelt2] * u.deg / u.pix)
    lofar_coord = SkyCoord(LOFAR_CORE.get_itrs(Time(obstime)))
    reference_coordinate = SkyCoord(header['crval1'] * u.deg, header['crval2'] * u.deg,
                                    frame='gcrs',
                                    observer=lofar_coord,
                                    obstime=obstime,
                                    distance=sun_distance)
    reference_coordinate_helio = reference_coordinate.transform_to(frames.Helioprojective)
    P1 = sun.P(obstime)
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
    area = (b_maj) * (b_min)*np.pi
    data *= area 
    data /= 1.e4
    image_map = sunpy.map.Map(data, header_heliocentric_frame)
    image_map = image_map.rotate()

    bl = SkyCoord(-fov * u.arcsec, -fov * u.arcsec, frame=image_map.coordinate_frame)
    tr = SkyCoord(fov * u.arcsec, fov * u.arcsec, frame=image_map.coordinate_frame)
    image_map = image_map.submap(bottom_left=bl, top_right=tr)

    return image_map

  
def image_plotter(fig: plt.Figure, ax: plt.Axes, sunpy_map: sunpy.map.Map):
    sunpy_map.draw_limb()
    solar_angle = sun.P(sunpy_map.date)
    scale = 3600 / 6000
    beam0 = mpl_patches.Ellipse((.7, -.7), sunpy_map.meta['bmaj'] * scale, sunpy_map.meta['bmin'] * scale,
                            angle=-(sunpy_map.meta['bpa']*u.deg + 90*u.deg).value, color='w', transform=ax.get_transform('world'))
    
    mapsmax = np.percentile(sunpy_map.data, q=99)
    mapsmin = sunpy_map.min()
    sunpy_map.plot_settings['cmap'] = plt.get_cmap('inferno')
    # sunpy_map.plot_settings['norm'] = mpl_colors.Normalize(vmin=mapsmin, vmax=mapsmax)

    sunpy_map.plot(axes=ax)
