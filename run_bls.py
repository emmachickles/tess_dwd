output_dir = '/scratch/echickle/bls/'
img_dir = output_dir + 'imgs/'
data_dir = '/scratch/data/tess/lcur/spoc/raws/'

# DWD: ./sector-05/tess2018319095959-s0005-0000000471013547-0125-s_lc.fits
# https://arxiv.org/pdf/2106.15104.pdf

# from __init__ import *
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import pdb
import os
import sys
import gc
sys.path.append(os.getcwd())

from KBB_Utils.KBB_Utils.Period_Finding import *

from astropy.io import fits
from astropy.timeseries import TimeSeries
from astropy.time import Time
from astropy.coordinates import SkyCoord
import astropy.units as u
import lc 
import pandas as pd 

# >> load white dwarf catalog
wd_cat  = pd.read_csv('/home/echickle/work/tess_dwd/WDs.txt', header=None, sep='\s+')

# >> WDs observed in S30 short cadence
# s30 = np.loadtxt('/home/echickle/work/tess_dwd/all_targets_S030_v1.txt')
# n = len(np.intersect1d(s30[:,0], wd_cat[0]))

# >> load flux and time for white dwarfs in S30 (from TESS_Photometry.py)
flux = np.load('/scratch/echickle/dwd/lc.npy').T
time = np.load('/scratch/echickle/dwd/ts.npy') 
ticid = np.load('/scratch/echickle/dwd/id.npy').astype('int')

# >> compute BLS
for i in range(len(flux)):
    y = flux[i]
    inds = np.nonzero(~np.isnan(y))
    t, y = time[inds], y[inds]
    dy = np.ones(y.shape)
    # default: pmin=3,pmax=True,qmin=2e-2,qmax=0.12
    t, y, dy, period, bls_power_best = BLS(t,y,dy,pmin=3,pmax=True,qmin=2e-2,qmax=0.3,remove=True)

    # >> compute phase curve
    prefix = str(bls_power_best)+'-'+str(round(period,2))+'-TIC%016d'%ticid[i]+'_'
    lc.make_phase_curve(t, y, period, output_dir=img_dir, prefix=prefix)

pdb.set_trace()

# fnames = []
# for i in os.listdir(data_dir):
#     for j in os.listdir(data_dir+i):
#         fnames.append(data_dir+i+'/'+j)

# with open(output_dir+'bls_power.txt', 'w') as f:
#     f.write('Sector,TICID,period,bls_power_best\n')

# # for fname in [data_dir + 'sector-05/tess2018319095959-s0005-0000000471013547-0125-s_lc.fits']:
# for fname in fnames:
#     try:
#         with fits.open(fname, memmap=False) as hdul:
#             y = hdul[1].data['PDCSAP_FLUX'] # >> unit = e-/s
#             dy = hdul[1].data['PDCSAP_FLUX_ERR'] # >> unit = e-/s
#             t = hdul[1].data['TIME'] # >> unit = BJD - 2457000, days

#     except:
#         print('Failed to open the following FITS file:')
#         print(fname)
#     gc.collect()

#     # >> mask out nans
#     inds = np.nonzero(~np.isnan(y))
#     y, dy, t = y[inds], dy[inds], t[inds]

#     # >> compute BLS
#     t, y, dy, period, bls_power_best = BLS(t,y,dy,pmin=3,pmax=True,qmin=2e-2,qmax=0.12,remove=True)
#     # >> pmin: minimum period in minutes, default 3 minutes
#     # >> pmax: maximum period in days, default 4/baseline

#     s = int(fname.split('-')[2][1:])
#     tic = int(fname.split('-')[3])
#     with open(output_dir+'bls_power.txt', 'a') as f:
#         f.write(str(s)+','+str(tic)+','+str(period)+','+str(bls_power_best)+'\n')


#     # >> compute phase curve
#     ts = TimeSeries(time=Time(t, format='jd'), data={'flux': y})
#     ts_folded = ts.fold(period=period*u.d) 

#     sorted_inds = np.argsort(ts_folded['time'])
#     folded_t = ts_folded['time'][sorted_inds]
#     folded_y = ts_folded['flux'][sorted_inds]

#     plt.figure()
#     plt.plot(folded_t.value, folded_y.value, '.k', ms=1)
#     plt.xlabel('Phase')
#     plt.ylabel('Flux')
#     plt.savefig(img_dir+str(bls_power_best)+'-'+fname.split('-')[2]+'-'+fname.split('-')[3]+'-'+str(period)+'_folded.png')
