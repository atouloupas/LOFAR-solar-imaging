import os
import numpy as np
import astropy.io.fits as fits
import astropy.wcs as wcs_lib
import matplotlib.pyplot as plt


def plot_fits_image(filename: str, output_path: str, with_colorbar: bool = True) -> str:
    """Plot a single FITS image in RA/Dec."""
    hdu = fits.open(filename)[0]
    wcs = wcs_lib.WCS(hdu.header)
    image = hdu.data[0, 0, :, :]

    fig = plt.figure(dpi=200)
    ax = plt.subplot(projection=wcs[0, 0, :, :])
    im = plt.imshow(image, cmap='inferno')
    plt.xlabel('Right ascension')
    plt.ylabel('Declination')

    if with_colorbar:
        fig.colorbar(im, ax=ax, pad=0.01)

    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    plt.savefig(output_path, bbox_inches='tight')
    plt.close(fig)
    return output_path
