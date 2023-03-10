from astropy.io import fits
import lc_utils as lcu
import pdb
import numpy as np
import matplotlib.pyplot as plt
from wotan import flatten
import sys
sys.path.insert(0, "/home/submit/echickle/work/")    
from KBB_Utils.KBB_Utils.Period_Finding import BLS

n_std=5
wind=0.1

# hdul = fits.open('/data/submit/echickle/foo/lc/hlsp_tglc_tess_ffi_gaiaid-1629388752470472704-s0014-cam3-ccd1_tess_v1_llc.fits')
hdul = fits.open('/data/submit/echickle/foo1/lc/hlsp_tglc_tess_ffi_gaiaid-1629388752470472704-s0041-cam3-ccd4_tess_v1_llc.fits')

t = hdul[1].data['time']
y = hdul[1].data['cal_psf_flux']
per = 49.7/1440
# per = 1/20.727432959364357

# >> remove nans
inds = np.nonzero(~np.isnan(y))
t, y = t[inds], y[inds]

# >> sigma-clip
med = np.median(y)
std = np.std(y)
inds = np.nonzero( (y > med - n_std*std) * (y < med + n_std*std) )
t, y = t[inds], y[inds]

# >> normalize        
if np.median(y) < 0:
    y = y / np.abs(med) + 2.
else:
    y = y / med

# >> detrending 
y = flatten(t, y, window_length=wind, method='biweight')
inds = np.nonzero(~np.isnan(y))
t, y = t[inds], y[inds]

# >> BLS
dy = np.ones(y.shape)
_,_,_, period, bls_power_best, freqs, power, dur, epo = \
    BLS(t,y,dy,pmin=20,pmax=0.25,qmin=0.005,qmax=0.2,remove=False)

# >> fold
t = t % per * 1440
inds = np.argsort(t)
t, y = t[inds], y[inds]

# >> bin
dy = np.ones(y.shape)
t, y, dy = lcu.bin_timeseries(t, y, 60, dy=dy)

# >> plot
shift = np.max(t) - np.min(t)
plt.figure(figsize=(5,3))
plt.errorbar(t, y, yerr=np.ones(len(y))*0.1,
               fmt='.k', ms=1, elinewidth=1)
plt.errorbar(t+shift, y, yerr=np.ones(len(y))*0.1,
               fmt='.k', ms=1, elinewidth=1)
plt.xlabel('Time [minutes]')
plt.ylabel('Relative Flux')
plt.tight_layout()
plt.savefig('/home/submit/echickle/gaia14aae_tglc_phase_curve.png', dpi=300)

plt.figure(figsize=(4,2))
#plt.plot(freqs, power, '-k', lw=0.5, alpha=0.7)
plt.plot(freqs, power, '.k', ms=0.5, alpha=0.7)
plt.xlabel('Frequency [1/days]')
plt.ylabel('BLS Power')
plt.title('period: '+str(round(period*1440,2))+' min')
plt.savefig('/home/submit/echickle/gaia14aae_tglc_spectrum.png')

plt.figure(figsize=(5,3))
plt.title('period :{}'.format(np.round(period*1440,2)))
plt.plot(1440/freqs, power, '-k', lw=0.5, alpha=0.7)
plt.xlim([0,120])
plt.xlabel('Period [minutes]')
plt.ylabel('BLS Power')
plt.tight_layout()
plt.savefig('/home/submit/echickle/gaia14aae_tglc_spectrum_zoom.png', dpi=300)

plt.figure(figsize=(4,2))
plt.hist(power, bins=150)
plt.xlabel('Power')
plt.ylabel('Frequency bins')
plt.savefig('/home/submit/echickle/gaia14aae_tglc_hist.png')

plt.figure(figsize=(4,2))
plt.hist(power,bins=500)
plt.xlabel('Power')
plt.ylabel('Frequency bins')
plt.xscale('log')
plt.savefig('/home/submit/echickle/gaia14aae_tglc_hist.png')


counts, bins, _ = plt.hist(power,bins=500)
p = bins[np.argmax(counts)]

# t, y, flag = lcu.prep_lc(t, y, n_std=3, wind=0.1)
# lcu.plot_phase_curve(t,y, per, '/home/submit/echickle/')
