gaia_tab = "/scratch/echickle/100pc_clean.fits"
out_dir = "/home/echickle/"

import pdb
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import astropy.units as u
from astropy.coordinates import SkyCoord
from astroquery.gaia import Gaia
import lc_utils as lcu

plt.rcParams['font.family'] = 'serif'
plt.rcParams["errorbar.capsize"] = 1
fig, ax = plt.subplots(figsize=(4,3))

# -- all gaia ------------------------------------------------------------------
    
hdul = fits.open(gaia_tab)
gid = hdul[1].data['source_id']
bp_rp = hdul[1].data['bp_rp']

parallax = hdul[1].data['parallax']
gmag = hdul[1].data['phot_g_mean_mag']
abs_mag = gmag+5*(np.log10(parallax)-2)

ax.plot(bp_rp, abs_mag, '.k', alpha=0.2, ms=0.05)
ax.set_xlim([-0.6, 5])
ax.set_ylim([17.5, 0])
ax.set_xlabel('Gaia BP-RP')
ax.set_ylabel('Absolute Magnitude (Gaia G)')

# -- CVs -----------------------------------------------------------------------

vmin = 60
vmax = 190
cm = plt.cm.get_cmap('plasma')

ticid = [800042858, 455206965,826164830,193092806] # >> low color
periods = [93.73902, 138.58632,108.65434,170.65144]
bp_rp = [-0.31, -0.39,-0.02,-0.37]
abs_mag = [8.27, 4.28,9.86,4.2]
sc = ax.scatter(bp_rp, abs_mag, marker="o", c=periods, cmap=cm, s=40,
                vmin=vmin, vmax=vmax, linewidths=0.7, edgecolors='k',
                alpha=1)

# CV candidates
ticid = [452954413,767706310,96852226,677736827,5393020,875850017,
         874726613,874383420,192991819,192991819,372519345,272551828]
periods = [113.94125, 94.16457,105.96957,98.07214,95.70348,85.12349,
           111.00473,86.21915,108.6033,108.60333,90.89413,107.27908]
bp_rp = [0.55, 0.44,0.75,0.64,0.44,0.33,
         0.49,0.38,0.59,0.59,0.43,0.75]
abs_mag = [10.24, 9.98,10.77,10.99,10.54,12.34,
           11.71,10.05,9.34,9.34,10.83, 10.54]
sc = ax.scatter(bp_rp, abs_mag, marker="o", c=periods, cmap=cm, s=40,
                vmin=vmin, vmax=vmax, linewidths=0.7, edgecolors='k',
                alpha=1)

ticid = [803489769, 36085812, 101433897, 808364853]
periods = [174.20556, 162.49021, 148.82092, 126.63071] # >> PERIOD GAP
bp_rp = [-0.05, 0.55, 0.69, 0.67]
abs_mag = [11.55, 12.1, 12.14, 9.61]
sc = ax.scatter(bp_rp, abs_mag, marker="v", c=periods, cmap=cm, s=40,
                vmin=vmin, vmax=vmax, linewidths=0.7, edgecolors='k',
                alpha=1)


cbar = fig.colorbar(sc, label="Orbital period (minutes)")

# -- save figure ---------------------------------------------------------------

plt.tight_layout()
# plt.savefig(out_dir+'hr.png', dpi=300)
plt.savefig(out_dir+'hr.pdf')
print(out_dir+'hr.png')

# -- phase curves --------------------------------------------------------------

sector = [61, 61, 62, 62]
cam = [1, 1, 1, 3]
ccd = [1, 1, 2, 4]
gid = [3070712053964938624, 3043990932113003648, 5670594640294856576, 5302515565071142400]
tess_periods = [174.2141, 162.49364, 148.82258, 126.63298]
bins = 40
lc_dir = "/scratch/echickle/grnd_lc/"

fig, ax = plt.subplots(nrows=2, ncols=4, figsize=(9,3))

for i in range(len(periods)):
    data_dir = "/scratch/data/tess/lcur/ffi/s%04d-lc/"%sector[i]
    suffix = '-'+str(cam[i])+'-'+str(ccd[i])+'.npy'
    ticid_list = np.load(data_dir+'id'+suffix).astype('int')
    ind = np.nonzero(ticid_list == ticid[i])[0][0]
    y = np.load(data_dir+'lc'+suffix)[ind]
    coord = np.load(data_dir+'co'+suffix)[ind]
    t = np.load(data_dir+'ts'+suffix)     

    y, _ = lcu.normalize_lc(y)    
    # t, y, flag = lcu.prep_lc(t, y, n_std=n_std, detrend=detrend, wind=wind)
    
    period = tess_periods[i]/1440
    folded_t, folded_y, folded_dy = lcu.bin_timeseries(t%period, y, bins)
    lcu.plot_phase_curve(ax[0][i], folded_t, folded_y, folded_dy, period,
                         ylabel=None, text_period=False)
    ax[0][i].text(0.95, 0.97, str(np.round(period*1440,5))+" min",
                  horizontalalignment="right", transform=ax[0][i].transAxes, va='top') 
    ax[0][i].set_xticklabels([])
    ax[0][i].set_xlabel('')
    lim = list(ax[0][i].get_ylim())
    lim[1] += 0.05
    ax[0][i].set_ylim(lim)

    t, y, dy, _, _ = lcu.load_atlas_lc(lc_dir + str(gid[i]))
    y, dy = lcu.normalize_lc(y, dy)
    period = periods[i]/1440
    folded_t, folded_y, folded_dy = lcu.bin_timeseries(t%period, y, bins, dy=dy)
    lcu.plot_phase_curve(ax[1][i], folded_t, folded_y, folded_dy, period,
                         ylabel=None, text_period=False)
    ax[1][i].text(0.95, 0.97, str(np.round(period*1440,5))+" min",
                  horizontalalignment="right", transform=ax[1][i].transAxes, va='top')
    lim = list(ax[1][i].get_ylim())
    lim[1] += 0.15
    ax[1][i].set_ylim(lim)


    # if i == 0:
    #     ax[1][i].set_ylim(ax[0][i].get_ylim())
    # elif i == 1:
    #     ax[0][i].set_ylim(ax[1][i].get_ylim())
    # elif i == 2:
    #     ax[0][i].set_ylim(ax[1][i].get_ylim())
    # elif i == 3:
    #     ax[0][i].set_ylim(ax[1][i].get_ylim())        


ax[0][0].set_ylabel('TESS\nNormalized Flux')
ax[1][0].set_ylabel('ATLAS\nNormalized Flux')


fig.tight_layout()
# plt.subplots_adjust(wspace=0.1, hspace=0.1)
# fig.savefig(out_dir+'phase_curve.png', dpi=300)
fig.savefig(out_dir+'phase_curve.pdf')
print(out_dir+'phase_curve.png')    




