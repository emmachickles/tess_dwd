# -- inputs --------------------------------------------------------------------

n_std = 3
wind  = 0.5
pmin = 400 / 60 # Nyquist in minutes
pmax = 0.25 # days

nbatches = 77
freq_batch_size = 14076

output_dir = '/home/echickle/out/'
data_dir   = '/home/echickle/data/s0061/s0061-lc/'

bls_dir    = output_dir + 's0061-bls-230228/'
ls_dir    = output_dir + 's0061-ls-230228/'
diag_dir   = output_dir + 's0061-diag-230228/'

# ------------------------------------------------------------------------------

import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt

from Period_Finding import frequency_grid, remove_harmonics

import pdb
import os
import sys
import gc
import fnmatch
from astropy.io import fits
from astropy.timeseries import BoxLeastSquares, LombScargle
from astropy.timeseries import TimeSeries
from astropy.time import Time
from astropy.coordinates import SkyCoord
import astropy.units as u
from wotan import flatten
import lc_utils as lcu
import gc

cam = sys.argv[1]
ccd = sys.argv[2]

os.makedirs(bls_dir, exist_ok=True)
os.makedirs(ls_dir, exist_ok=True)
os.makedirs(diag_dir, exist_ok=True)

fail_txt = output_dir + 'cam'+str(cam)+'-failed.txt'
with open(fail_txt, 'w') as f:
    f.write('')

suffix = '-'+str(cam)+'-'+str(ccd)+'.npy'
flux = np.load(data_dir+'lc'+suffix)
time = np.load(data_dir+'ts'+suffix)     
coord = np.load(data_dir+'co'+suffix)
ticid = np.load(data_dir+'id'+suffix).astype('int')

inds = np.argsort(time)
time, flux = time[inds], flux[:,inds]

# >> compute BLS
for i in range(len(flux)):
    print(i)
    y = flux[i]
    t = time

    # -- prep light curve --------------------------------------------------
    t, y, flag = lcu.prep_lc(t, y, n_std=n_std, wind=wind, diag=False,
                             ticid=ticid[i], cam=cam, ccd=ccd,
                             coord=coord[i], output_dir=diag_dir)
    dy = np.ones(y.shape)
    
    if flag:
        with open(fail_txt, 'a') as f:
            f.write(str(ticid[i])+'\n')
    else:

        # -- compute BLS ---------------------------------------------------
        freqs = frequency_grid(t,y,pmin=pmin,pmax=pmax, qmin=0.005, qmax=0.05)
        periods = 1/freqs * u.day
        model = BoxLeastSquares(t * u.day, y, dy)

        # # >> code inspired by cuvarbase
        for batch in range(nbatches):
            imin = freq_batch_size * batch
            imax = min([len(freqs), freq_batch_size * (batch + 1)])
            minp = 1/freqs[imax] # >> min period,
        # maximum transit duration must be shorter than min period
        minq = locext(min, qmin, imin, imax)
        maxq = locext(max, qmax, imin, imax)            
            
        periodogram = model.power(periods, 0.005)
        pdb.set_trace()
        freqs, power, dur, epo = remove_harmonics(freqs, periodogram.power, periodogram.duration,
                                        periodogram.transit_time)
        max_power = np.argmax(power)
        period, bls_power_best = 1/freqs[max_power], power[max_power]
        dur, epo = dur[max_power], epo[max_power]
        pdb.set_trace()
        
        # -- compute LS ----------------------------------------------------
        ls_power = LombScargle(t, y).power(freq)
        max_power = np.argmax(ls_power)
        ls_power, ls_power_best = 1/freqs[max_power], ls_power[max_power]
        
        # -- plot phase curve ----------------------------------------------
        prefix = 'pow_'+str(bls_power_best)+'_per_'+str(round(period*1440,5))+\
            '_TIC%016d'%ticid[i]+'_cam_'+str(cam)+'_ccd_'+str(ccd)+\
            '_dur_'+str(dur)+'_epo_'+str(epo)+\
            '_ra_{}_dec_{}_'.format(coord[i][0], coord[i][1])                

        lcu.make_phase_curve(t, y, period, dy=dy, output_dir=bls_dir,
                             prefix=prefix, freqs=freqs, power=power,
                             ticid=ticid[i], bins=100)

        prefix = 'pow_'+str(ls_power_best)+'_per_'+str(round(ls_period*1440,5))+\
            '_TIC%016d'%ticid[i]+'_cam_'+str(cam)+'_ccd_'+str(ccd)+\
            '_ra_{}_dec_{}_'.format(coord[i][0], coord[i][1])                

        lcu.make_phase_curve(t, y, ls_period, dy=dy, output_dir=ls_dir,
                             prefix=prefix, freqs=freqs, power=ls_power,
                             ticid=ticid[i], bins=100, bls=False)
        

