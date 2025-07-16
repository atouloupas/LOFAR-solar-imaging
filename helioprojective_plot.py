from astropy.io import fits
import matplotlib.pyplot as plt
import matplotlib.colors as mpl_colors
import matplotlib.patches as mpl_patches
import numpy as np
import astropy.units as u
from imaging_utils import reproject_to_heliocentric_frame
import os
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("filename", help="FITS file")
args = parser.parse_args()

# read in fits file and pull out header and data
file = args.filename
basename = os.path.splitext(os.path.basename(file))[0]
path_out = f"images_png/helioprojective_{basename}.png"

hdu = fits.open(file)
header = hdu[0].header
data = np.squeeze(hdu[0].data)  # squeezed here as data in shape [1,1,1500,1500]

fov = 3000
rotate_map = reproject_to_heliocentric_frame(header, data, fov)
fig = plt.figure(dpi=200)
ax = plt.subplot(projection=rotate_map)

rotate_map.draw_limb()
rotate_map.plot_settings['cmap'] = plt.get_cmap('inferno')
mapsmin = 0.
mapsmax = rotate_map.max()
rotate_map.plot_settings['norm'] = mpl_colors.Normalize(vmin=mapsmin, vmax=mapsmax)
scale = 3600 / 6000
beam0 = mpl_patches.Ellipse((.7, -.7), rotate_map.meta['bmaj'] * scale, rotate_map.meta['bmin'] * scale,
                            angle=-(rotate_map.meta['bpa']*u.deg + 90*u.deg).value, color='w', transform=ax.get_transform('world'))
ax.add_patch(beam0)

im = rotate_map.plot(axes=ax)
cbar = fig.colorbar(im, ax=ax, pad=0.01)
cbar.ax.text(0.4, 0.5, "Flux [SFU/pixel]", rotation=270, color='w', fontsize='medium',
             ha='center', va='center', transform=cbar.ax.transAxes)
fig.savefig(path_out, bbox_inches='tight')
plt.close()

print(f"Helioprojective image completed: {path_out}")