from casacore.tables import table, taql
import scipy as sp
import numpy as np
from astropy.time import Time
from pathlib import Path
import argparse
import matplotlib.pyplot as plt
from matplotlib.legend_handler import HandlerTuple

parser = argparse.ArgumentParser()
parser.add_argument("filename", help="MS file")
parser.add_argument("--std_thresh", default=2, type=float, help="Standard deviation threshold for outliers")
parser.add_argument("--remote", action="store_true", help="Plot remote stations")
parser.add_argument("--plot", action="store_true", help="Plot file")
args = parser.parse_args()

file = args.filename
filename = Path(file)
while filename.suffix:
    filename = filename.with_suffix('')
filename = filename.name

with table(file) as dataset:
    t = taql("select ANTENNA1, DATA, TIME from $dataset where ANTENNA1=ANTENNA2")
    antenna1 = t.getcol("ANTENNA1")
    antenna2 = antenna1
    time = t.getcol("TIME")
    visibilities_ac = t.getcol("DATA")

with table(file + "/SPECTRAL_WINDOW") as spectral_window:
    frequencies = spectral_window[0]["CHAN_FREQ"]  # channel/subband frequencies
    
with table(file + "/ANTENNA") as antenna:
    ant = antenna.getcol("NAME")
    ant_idx = antenna2

if args.remote == False:
    ant = [name for name in ant if name.startswith("CS")]  # keep only core stations
    ant_idx = [idx for idx, name in enumerate(ant) if name.startswith("CS")]

# Get only values for selected antennas
mask = np.isin(antenna1, ant_idx)
visibilities_ac = visibilities_ac[mask]
time = time[mask]

num_chan = len(frequencies)  # number of frequency channels
sb_num = num_chan // 2  # subband number
frequency = frequencies[sb_num]  # frequency of specific subband

intensity_ac = np.sqrt(np.abs(visibilities_ac)) * 1e-4  # convert from Jy to SFU


# Convert MJD to UTC
def mjd_to_utc(time, str_format):
    """Convert MJD time (in seconds) to UTC.

    Parameters
    ----------
    time : array, float
        Array or value of time in MJD.
    str_format : str
        String format for displaying UTC.

    Returns
    -------
    date_time : array, str
        The times in the specified string format.
    """
    time_mjd = time / (24 * 3600)
    time_array = Time(time_mjd, format='mjd')
    date_time = time_array.strftime(str_format)

    return date_time

# Find the start time of the observation
date_time = mjd_to_utc(time[0], '%H:%M:%S')
start_time = mjd_to_utc(time[0], "%Y-%m-%d")

# Find intensity average for each time window
unique_times_ac = mjd_to_utc(np.unique(time), '%H:%M:%S')  # if not using taql, use mask
intensity_ac_avg = np.mean(intensity_ac[:,sb_num,0].reshape(len(unique_times_ac), -1), axis=1)

# Average intensity of all antennas
ant_int_avg = np.mean(intensity_ac_avg)
# Average intensity of each antenna
ant_int = [np.mean(intensity_ac[antenna1[mask] == i][:,sb_num,0]) for i in range(len(ant))]
# Standard deviations of the intensities
ant_int_std = np.std(ant_int)

# Thresholds
std_thresh = args.std_thresh  # number of standard deviations
upper_thresh = ant_int_avg + std_thresh * ant_int_std
lower_thresh = ant_int_avg - std_thresh * ant_int_std
# Outliers
outliers = (ant_int > upper_thresh) | (ant_int < lower_thresh)
outliers_idx = np.where(outliers)[0]

# Get the names of the outlier antennas and save to file
outlier_stations = np.array(ant)[outliers]
print(f'Outlier stations: {outlier_stations}')

with open(f'stations/stations_{filename}.txt', 'w') as f:
    for station in outlier_stations:
        f.write(f"{station}\n")


# If specified, plot intensities of stations with outliers
if args.plot == True:
    fig, ax = plt.subplots(dpi=200)

    for i, name in enumerate(ant):
        # Mark the outlier stations
        if outliers[i] == True:
            c = 'r'
        else:
            c = 'k'

        # Mark the core and remote stations
        if name.startswith('CS'):
            marker = 'o'
        else:
            marker = '^'

        ax.scatter(i, ant_int[i], c=c, marker=marker)

    p1 = ax.axhline(y=ant_int_avg, color='r', linestyle='-')
    p2 = ax.axhline(y=upper_thresh, color='g', linestyle='--')
    ax.axhline(y=lower_thresh, color='g', linestyle='--')
    ax.axhspan(lower_thresh, upper_thresh, facecolor='g', alpha=0.1)
    p3 = ax.scatter([], [], c='k', marker='o')
    p5 = ax.scatter([], [], c='r', marker='o')
    if args.remote:
        p4 = ax.scatter([], [], c='k', marker='^')

    ax.set_title(f"{frequency*1e-6:.2f} MHz, {start_time} {date_time}")
    ax.set_ylabel("Station flux density [SFU]")
    ax.set_xlabel("Station name")
    ax.set_xticks(range(len(ant)), ant, rotation=90, fontsize='small')
    # Mark outlier antenna labels
    for i in outliers_idx:
        ax.get_xticklabels()[i].set_color('r')

    if args.remote:
        handles = [p1, p2, p3, p4, p5]
        labels = ['flux density\naverage', fr'{std_thresh}$\sigma$ threshold', 'core stations', 'remote stations', 'flagged stations']
    else:
        handles = [p1, p2, p3, p5]
        labels = ['flux density\naverage', fr'{std_thresh}$\sigma$ threshold', 'core stations', 'flagged stations']

    plt.tight_layout()
    plt.legend(handles, labels, handler_map={tuple: HandlerTuple(ndivide=None)}, fontsize="small")
    plt.savefig(f"./plots/ant_intensities_{filename}.png")
    plt.close()