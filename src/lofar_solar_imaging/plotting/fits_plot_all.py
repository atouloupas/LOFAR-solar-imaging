import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import astropy.io.fits as fits
import astropy.wcs as wcs_lib
import os
import dateutil
import argparse
import glob
import re

parser = argparse.ArgumentParser()
parser.add_argument("filename", help="FITS file or directory to plot")
parser.add_argument("--no_cbar", action="store_true", help="Do not show colourbar")
parser.add_argument("--plot_all", action="store_true", help="Plot all generated images")
args = parser.parse_args()

def load_image_data(fits_path, plot_all):
    hdu = fits.open(fits_path)[0]
    header = hdu.header
    wcs = wcs_lib.WCS(header)
    image = hdu.data[0, 0, :, :]
    iso_date = header["DATE-OBS"]
    obs_date = dateutil.parser.isoparse(iso_date).strftime("%Y-%m-%d %H:%M:%S")
    freq_min = header["CRVAL3"] / 1e6
    d_freq = header["CDELT3"] / 1e6
    freq_max = freq_min + d_freq
    if not plot_all:
        b_maj = header['BMAJ']
        b_min = header['BMIN']
        area = b_maj * b_min * np.pi
        image *= area
        image /= 1.e4  # Convert to SFU/pixel
        
    return image, wcs, header, freq_min, obs_date

if args.plot_all:
    directory = args.filename.rstrip('/')
    base = os.path.basename(directory)
    pattern = os.path.join(directory, f"{base}-*.fits")
    all_files = sorted(glob.glob(pattern))

    # Exclude files containing 't' followed by 4 digits
    fits_files = [f for f in all_files if not re.search(r't\d{4}', os.path.basename(f))]

    if len(fits_files) == 0:
        raise FileNotFoundError(f"No valid FITS files found in: {directory}")

    types = [os.path.splitext(os.path.basename(f))[0].split(f"{base}-")[1] for f in fits_files]
    types = [i.upper() if i == 'psf' else i.title() for i in types]

    layout = "AABB;CCDD;.EE."

    fig, axs = plt.subplot_mosaic(layout, figsize=(9, 12), dpi=200, subplot_kw={'projection': None},
                                  gridspec_kw={'wspace': -0.1, 'hspace': 0.3})
    label_map = ['A', 'B', 'C', 'D', 'E']

    for i, (fits_path, img_type) in enumerate(zip(fits_files, types)):
        image, wcs, header, freq_min, obs_date = load_image_data(fits_path, args.plot_all)

        label = label_map[i]
        ax = axs[label]

        # Replace plain Axes with WCS-aware one
        fig.delaxes(ax)
        axs[label] = fig.add_subplot(ax.get_subplotspec(), projection=wcs[0, 0, :, :])

        im = axs[label].imshow(image, cmap='inferno', origin='lower')
        axs[label].set_title(f"{img_type}", fontsize='x-large')
        axs[label].set_xlabel("Right ascension", fontsize='small')
        axs[label].set_ylabel("Declination", fontsize='small')
        axs[label].tick_params(axis='both', labelsize='small')
        if label in ['B', 'D']:
            axs[label].tick_params(axis='y', labelleft=False)
            axs[label].set_ylabel("")
        # axs[label].axis('off')

        if not args.no_cbar:
            cbar = fig.colorbar(im, ax=axs[label], fraction=0.047, pad=0.01)
            cbar.ax.tick_params(labelsize='small')
    # plt.suptitle(f"{freq_min:.2f} MHz {obs_date}", fontsize='x-large')
    plt.savefig(f"./images_png/{base}_all.png", bbox_inches="tight")
    plt.close()


else:
    filename = args.filename
    basename = os.path.splitext(os.path.basename(filename))[0]
    target = "Sun" if ("sun" in filename.lower() or "sap000" in filename.lower()) else "Cas A"

    if "model" in filename.lower():
        status = "Model"
    elif "calibrated" in filename.lower():
        status = "Calibrated"
    else:
        status = "Original"

    image, wcs, header, freq_min, obs_date = load_image_data(filename, True)

    fig = plt.figure(dpi=200)
    ax = plt.subplot(projection=wcs[0, 0, :, :])
    im = plt.imshow(image, cmap='inferno')
    plt.xlabel("Right ascension")
    plt.ylabel("Declination")
    plt.title(f"{status}", fontsize='x-large')

    if not args.no_cbar:
        cbar = fig.colorbar(im, ax=ax, pad=0.01)
        # cbar.ax.text(0.4, 0.5, "Flux [SFU/pixel]", rotation=270, color='w',
        #              fontsize='medium', ha='center', va='center', transform=cbar.ax.transAxes)

    plt.savefig(f"./images_png/{basename}.png", bbox_inches="tight")
    plt.close()
