from astropy.io import fits
import lc_utils as lcu
import pdb
import numpy as np
import matplotlib.pyplot as plt
from wotan import flatten
import sys
sys.path.insert(0, "/home/submit/echickle/work/")    
from KBB_Utils.KBB_Utils.Period_Finding import BLS

n_std=3
wind=0.05

data_dir = '/data/submit/echickle/data/sector_41_lc/'
suffix = '-3-4'
ticid = 1201247611
ticid_list = np.load(data_dir+'id'+suffix+'.npy').astype('int')
ind = np.nonzero(ticid_list == ticid)[0][0]
t = np.load('/data/submit/echickle/data/sector_41_lc/ts'+suffix+'.npy')
y = np.load('/data/submit/echickle/data/sector_41_lc/lc'+suffix+'.npy')[ind]


per = 49.7/1440

# try sector 41
# try tutorial phase folding

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
    BLS(t,y,dy,pmin=7,pmax=0.25,qmin=0.005,qmax=0.2,remove=False)

# >> fold
t = t % per * 1440
inds = np.argsort(t)
t, y = t[inds], y[inds]

# >> bin
dy = np.ones(y.shape)
t, y, dy = lcu.bin_timeseries(t, y, dy, 200)

# >> plot
plt.figure()
plt.plot(t, y, '.k', ms=1)
shift = np.max(t) - np.min(t)        
plt.plot(t+shift, y, '.k', ms=1) 
plt.xlabel('Time [minutes]')
plt.ylabel('Relative Flux')
plt.savefig('/home/submit/echickle/gaia14aae_aplc_phase_curve.png')

plt.figure()
plt.plot(freqs, power, '-k', lw=0.5, alpha=0.7)
plt.xlabel('Frequency [1/days]')
plt.ylabel('BLS Power')
plt.title('period: '+str(round(period*1440,2))+' min')
plt.savefig('/home/submit/echickle/gaia14aae_aplc_spectrum.png')

plt.figure()
plt.plot(1440/freqs, power, '-k', lw=0.5, alpha=0.7)
plt.xlim([0,120])
plt.xlabel('Period [minutes]')
plt.ylabel('BLS Power')
plt.savefig('/home/submit/echickle/gaia14aae_aplc_spectrum_zoom.png')


# t, y, flag = lcu.prep_lc(t, y, n_std=3, wind=0.1)
# lcu.plot_phase_curve(t,y, per, '/home/submit/echickle/')
