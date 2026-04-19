import numpy as np
import astropy.units as u
from astropy.coordinates import EarthLocation, SkyCoord
from astropy.time import Time
import sunpy.map
from sunpy.coordinates import frames, sun

LOFAR_CORE = EarthLocation.from_geocentric(
    x=3826560.737625 * u.m,
    y=461033.24679167 * u.m,
    z=5064904.97320833 * u.m,
)


def reproject_to_heliocentric_frame(header, data, fov=3000):
    """Reproject an equatorial FITS image into helioprojective coordinates."""
    obstime = header['date-obs']
    freq = header['crval3'] * u.Hz
    cdelt1, cdelt2 = abs(header['cdelt1']), abs(header['cdelt2'])
    sun_distance = sun.earth_distance(obstime)
    reference_pixel = u.Quantity([header['crpix1'], header['crpix2']] * u.pixel)
    data = data.squeeze()
    scale = u.Quantity([cdelt1, cdelt2] * u.deg / u.pix)
    lofar_coord = SkyCoord(LOFAR_CORE.get_itrs(Time(obstime)))
    reference_coordinate = SkyCoord(
        header['crval1'] * u.deg,
        header['crval2'] * u.deg,
        frame='gcrs',
        observer=lofar_coord,
        obstime=obstime,
        distance=sun_distance,
    )
    reference_coordinate_helio = reference_coordinate.transform_to(frames.Helioprojective)
    p_angle = sun.P(obstime)
    header_heliocentric_frame = sunpy.map.make_fitswcs_header(
        data,
        reference_coordinate_helio,
        reference_pixel=reference_pixel,
        scale=scale,
        rotation_angle=-p_angle,
        wavelength=freq.to(u.MHz).round(2),
        observatory='LOFAR',
    )

    for key in ('HISTORY', 'BMAJ', 'BMIN', 'BPA'):
        if key in header:
            header_heliocentric_frame[key] = header[key]

    area = header['BMAJ'] * header['BMIN'] * np.pi
    data = data * area / 1.0e4
    image_map = sunpy.map.Map(data, header_heliocentric_frame).rotate()

    bl = SkyCoord(-fov * u.arcsec, -fov * u.arcsec, frame=image_map.coordinate_frame)
    tr = SkyCoord(fov * u.arcsec, fov * u.arcsec, frame=image_map.coordinate_frame)
    return image_map.submap(bottom_left=bl, top_right=tr)
