from casacore.tables import table, taql
import scipy as sp
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoLocator, LogLocator
import matplotlib.dates as mdates
from matplotlib.legend_handler import HandlerTuple
from astropy.time import Time
from imaging_utils import model_flux
from pathlib import Path


file = "HBA_joined_cas_20230311.MS"
filename = Path(file)
while filename.suffix:
    filename = filename.with_suffix('')
filename = filename.name

# Retrieve target / calibrator name for plot titles
target = "Sun" if ("sun" or "sap000") in filename.lower() else "Cas A"


with table(file) as dataset:
    t = taql("select ANTENNA1, DATA, TIME from $dataset where ANTENNA1=ANTENNA2")
    # antenna1 = dataset.getcol("ANTENNA1")
    # antenna2 = dataset.getcol("ANTENNA2")
    # uvw = dataset.getcol("UVW")
    # time = dataset.getcol("TIME")

    # Assuming single channel, single polarization
    # visibilities = dataset.getcol("DATA")
    # visibilities_ac = visibilities[antenna1 == antenna2]  # autocorrelations

    antenna1 = t.getcol("ANTENNA1")
    antenna2 = antenna1
    time = t.getcol("TIME")
    visibilities_ac = t.getcol("DATA")

with table(file + "/SPECTRAL_WINDOW") as spectral_window:
    # print(spectral_window[0]['NAME'])
    frequencies = spectral_window[0]["CHAN_FREQ"]  # channel/subband frequencies
    
with table(file + "/ANTENNA") as antenna:
    ant = np.array(antenna.getcol("NAME"))
    print(antenna.nrows())
    print(ant)


with open(f'stations/stations_{filename}.txt', 'r') as f:
    out_stations = f.read().splitlines()
    ant_out_idx = [np.where(ant == name)[0] for name in out_stations]
    ant_idx = [idx for idx, name in enumerate(ant) if (name.startswith("CS") and name not in out_stations)]


# Choose which antennas to keep
# ant_idx = range(len(ant))  # keep all
# ant_idx = [idx for idx, name in enumerate(ant) if name.startswith("CS")]  # keep only core stations
# ant_idx = ant.index("CS011LBA")

# Get only values for selected antennas
mask = np.isin(antenna1, ant_idx)
visibilities_ac = visibilities_ac[mask]
time = time[mask]

num_chan = len(frequencies)  # number of frequency channels
sb_num = num_chan // 2  # subband number
frequency = frequencies[sb_num]  # frequency of specific subband
wavelength = (sp.constants.c) / frequency
# uvw /= wavelength
# pol = ["$XX$", "$XY$", "$YX$", "$YY$"]

# intensity = np.abs(visibilities)  # visibility amplitude
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


date_time = mjd_to_utc(time, '%H:%M:%S')
start_time = mjd_to_utc(time[0], "%Y-%m-%d")
date_time_ac = date_time#[antenna1 == antenna2]

# Find intensity average for each time window
unique_times_ac = mjd_to_utc(np.unique(time), '%H:%M:%S')  # if not using taql, use mask
intensity_ac_avg = np.mean(intensity_ac[:,sb_num,0].reshape(len(unique_times_ac), -1), axis=1)
# intensity_ac_avg = intensity_ac[:,sb_num,0].reshape(len(unique_times_ac), -1)[:,0,:]

print(len(unique_times_ac), len(intensity_ac_avg), intensity_ac.shape)

# Plot amplitudes for both polarisations
fig, ax = plt.subplots(dpi=200)

# for i in range(len(ant)):
#     ax.scatter(date_time_ac[antenna1 == i], intensity_ac[antenna1 == i][:,sb_num,0], s=0.1, alpha=0.7)
ax.scatter(date_time_ac, intensity_ac[:,sb_num,0], label="autocorrelated values", s=0.5, alpha=1)
ax.plot(unique_times_ac, intensity_ac_avg, label="avg autocorrelated values", color="r", marker='o', alpha=1)
ax.set_title(f"Flux density of {target} at {frequency*1e-6:.2f} MHz ({start_time})")
ax.set_xlabel("Time [UTC]")
ax.set_ylabel("Flux density [SFU]")
ax.xaxis.set_major_locator(AutoLocator())
fig.autofmt_xdate(ha="center", rotation=20)

plt.legend()
plt.tight_layout()
plt.savefig(f"./plots/intensities_{filename}.png")
plt.close()


##############################################################################################

# Average intensity of all antennas
ant_int_avg = np.mean(intensity_ac_avg)
# Average intensity of each antenna
ant_int = [np.mean(intensity_ac[antenna1[mask] == i][:,sb_num,0]) for i in ant_idx]
# Standard deviations of the intensities
ant_int_std = np.std(ant_int)

# Thresholds
upper_thresh = ant_int_avg + 2 * ant_int_std
lower_thresh = ant_int_avg - 2 * ant_int_std
# Outliers
outliers = (ant_int > upper_thresh) | (ant_int < lower_thresh)
outliers_idx = np.where(outliers)[0]

# Get the names of the outlier antennas and save to file
outlier_stations = ant[ant_idx][outliers]
with open(f'stations/stations_{filename}.txt', 'w') as f:
    for station in outlier_stations:
        f.write(f"{station}\n")

fig, ax = plt.subplots(dpi=200)

for i, name in enumerate(ant[ant_idx]):
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
p4 = ax.scatter([], [], c='k', marker='^')
p5 = ax.scatter([], [], c='r', marker='o')
p6 = ax.scatter([], [], c='r', marker='^')

ax.set_title(f"Station flagging ({frequency*1e-6:.2f} MHz {start_time} {date_time[0]})")
ax.set_ylabel("Station flux density [SFU]")
ax.set_xlabel("Station name")
ax.set_xticks(range(len(ant[ant_idx])), ant[ant_idx], rotation=90, fontsize='small')
# Mark outlier antenna labels
for i in outliers_idx:
    ax.get_xticklabels()[i].set_color('r')

handles = [p1, p2, p3, p4, (p5, p6)]
labels = ['flux density\naverage', r'2$\sigma$ threshold', 'core stations', 'remote stations', 'outlier stations']

plt.tight_layout()
plt.legend(handles, labels, handler_map={tuple: HandlerTuple(ndivide=None)}, fontsize="small")
plt.savefig(f"./plots/ant_intensities_{filename}.png")
plt.close()

##############################################################################################

# Convert MJD to UTC and get mdates for plotting
utc_times = Time(np.array([time[0], time[-1]]) / (24*3600), format='mjd').to_datetime()
date_nums = mdates.date2num(utc_times)

# COMBINE FREQUENCIES
intensity_ac_avg_all = np.mean(intensity_ac[:,:,0].reshape(len(unique_times_ac), -1, num_chan), axis=1)

# Get average intensity in time for each channel and subtract from data
transmission = np.nanmean(intensity_ac_avg_all, axis=0)  # this is T(nu)

# Get model flux of Cas A for correction
flux_model = model_flux('CasA', frequencies)  # this is M_CasA(nu)
bandpass = transmission / flux_model  # this is B(nu)

# Get sky flux: S_sky(t,nu) = F(t,nu) / B(nu) - M_CasA(nu)
sky_flux = intensity_ac_avg_all / bandpass - flux_model

# print(chan_intensity)
fig, ax = plt.subplots(dpi=200)
ax.plot(frequencies*1e-6, transmission, label='average total flux',
         linewidth=1, marker='o', markersize=2, alpha=0.5)
ax.plot(frequencies*1e-6, bandpass, label='bandpass',
         linewidth=1, marker='o', markersize=2, alpha=0.5)
ax.plot(frequencies*1e-6, flux_model, label='Cas A model flux')
ax.set_xlabel('Frequency [MHz]')
ax.set_ylabel('Flux density [SFU]')
# ax.set_title(f"Bandpass of {target}")
ax.set_yscale('log')
plt.legend()
plt.tight_layout()
plt.savefig(f"./plots/bandpass_{filename}.png")
plt.close()


###############################################################################################

# Plot dynamic spectrum
percentiles = np.nanpercentile(sky_flux, [20, 90])
vmin, vmax = percentiles
# vmin = max(vmin, 0.1)
# vmax = min(vmax, 7)

fig, ax = plt.subplots(dpi=200)

img = ax.imshow(sky_flux.T, aspect='auto', cmap='inferno', 
                norm='linear', origin='lower',
                # vmin=vmin, vmax=vmax,
                vmin=-.02, vmax=0.03,
                extent=[date_nums[0], date_nums[-1],
                        frequencies[0]*1e-6, frequencies[-1]*1e-6])
ax.set_title(f"Dynamic spectrum of {target} ({start_time})")
ax.set_xlabel("Time [UTC]")
ax.set_ylabel("Frequency [MHz]")
ax.xaxis_date()
ax.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M:%S"))
plt.gcf().autofmt_xdate(rotation=20, ha='center')
# ax.xaxis.set_major_locator(mdates.MinuteLocator(byminute=[0, 15, 30, 45], interval=1))
plt.minorticks_on()
ax.xaxis.set_major_locator(AutoLocator())

cbar = fig.colorbar(img, pad=0.01)
cbar.ax.text(0.4, 0.5, "Flux density [SFU]", rotation=270, color='w', fontsize='medium',
             ha='center', va='center', transform=cbar.ax.transAxes)

plt.tight_layout()
plt.savefig(f"./plots/dynspec_{filename}.png")
plt.close()