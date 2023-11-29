import matplotlib as mpl
import matplotlib.pyplot as plt
from astropy.io import fits
import numpy as np
import pandas as pd
from matplotlib.colors import LogNorm

wd_tab = "/home/echickle/data/WDs.txt"
wd_main = "/data/GaiaEDR3_WD_main.fits"
rp_ext = "/data/GaiaEDR3_WD_RPM_ext.fits"

source_id = np.empty(0)
bp_rp = np.empty(0)
parallax = np.empty(0)
gmag = np.empty(0)

maincat = fits.open(wd_main)
source_id = np.append(source_id, maincat[1].data['source_id'])
bp_rp = np.append(bp_rp, maincat[1].data['bp_rp'])
parallax = np.append(parallax, maincat[1].data['parallax'])
gmag = np.append(gmag, maincat[1].data['phot_g_mean_mag'])

rpmext = fits.open(rp_ext) 
source_id = np.append(source_id, rpmext[1].data['source_id'])
bp_rp = np.append(bp_rp, rpmext[1].data['bp_rp'])
parallax = np.append(parallax, rpmext[1].data['parallax'])
gmag = np.append(gmag, rpmext[1].data['phot_g_mean_mag'])

abs_mag = gmag+5*(np.log10(parallax)-2)


# Set the default text color to white
mpl.rcParams['text.color'] = 'white'
mpl.rcParams['axes.labelcolor'] = 'white'
mpl.rcParams['xtick.color'] = 'white'
mpl.rcParams['ytick.color'] = 'white'

# Create a figure and axis with a black background
fig, ax = plt.subplots(figsize=(5, 4))
fig.patch.set_facecolor('black')
ax.set_facecolor('black')

hist = ax.hist2d(bp_rp, abs_mag, bins=200, range=[[-0.6, 1.9], [4, 15.5]], density=False,
                 cmin=0.03, norm=LogNorm(), cmap='plasma')
# ax.tick_params(top=True, labeltop=True, bottom=False, labelbottom=False)
ax.set_xlabel('Gaia BP-RP', color='white')
ax.set_ylabel('Absolute Magnitude (Gaia G)', color='white')
ax.invert_yaxis()

# Add a colorbar
cbar = plt.colorbar(hist[3], ax=ax)
cbar.ax.set_ylabel('Counts', color='white')
cbar.ax.yaxis.set_tick_params(color='white')
cbar.ax.xaxis.set_tick_params(color='white')
plt.setp(plt.getp(cbar.ax.axes, 'yticklabels'), color='white')

plt.savefig('/home/echickle/hr.png', dpi=300)
print('/home/echickle/hr.png')
