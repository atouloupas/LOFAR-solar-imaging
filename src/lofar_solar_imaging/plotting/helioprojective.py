import os
import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
import matplotlib.colors as mpl_colors
import matplotlib.patches as mpl_patches
import astropy.units as u

from lofar_solar_imaging.reprojection import reproject_to_heliocentric_frame


def plot_helioprojective_fits(filename: str, output_path: str, fov: int = 3000) -> str:
    """Generate a helioprojective PNG from one FITS image."""
    with fits.open(filename) as hdul:
        header = hdul[0].header
        data = np.squeeze(hdul[0].data)

    rotate_map = reproject_to_heliocentric_frame(header, data, fov=fov)
    fig = plt.figure(dpi=200)
    ax = plt.subplot(projection=rotate_map)

    rotate_map.draw_limb()
    rotate_map.plot_settings['cmap'] = plt.get_cmap('inferno')
    rotate_map.plot_settings['norm'] = mpl_colors.Normalize(vmin=0.0, vmax=rotate_map.max())
    scale = 3600 / 6000
    beam0 = mpl_patches.Ellipse(
        (.7, -.7),
        rotate_map.meta['bmaj'] * scale,
        rotate_map.meta['bmin'] * scale,
        angle=-(rotate_map.meta['bpa'] * u.deg + 90 * u.deg).value,
        color='w',
        transform=ax.get_transform('world'),
    )
    ax.add_patch(beam0)

    image = rotate_map.plot(axes=ax)
    cbar = fig.colorbar(image, ax=ax, pad=0.01)
    cbar.ax.text(0.4, 0.5, 'Flux [SFU/pixel]', rotation=270, color='w', fontsize='medium',
                 ha='center', va='center', transform=cbar.ax.transAxes)

    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    fig.savefig(output_path, bbox_inches='tight')
    plt.close(fig)
    return output_path
