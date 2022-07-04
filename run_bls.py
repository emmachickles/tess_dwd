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
# import lc 
import pandas as pd 


# lcfile = '/scratch/submit/tess/data/tesscurl_sector_5_lc/tess2018319095959-s0005-0000000471013547-0125-s_lc.fits'

lcfile = '/scratch/echickle/tmp/tess2021204101404-s0041-0000001201247611-0212-s_lc.fits'
hdul = fits.open(lcfile)
t = hdul[1].data['TIME']
y = hdul[1].data['PDCSAP_FLUX']

# >> remove quality flags
inds = np.nonzero(hdul[1].data['QUALITY'])
t = np.delete(t, inds)
y = np.delete(y, inds)

inds = np.nonzero(~np.isnan(y))
t, y = t[inds], y[inds]
dy = np.ones(y.shape)
plt.figure()
plt.plot(t, y, '.k', ms=1)
plt.xlabel('TIME')
plt.ylabel('PDCSAP_FLUX')
plt.title('TIC 1201247611, S41-20s')
plt.savefig('/scratch/echickle/tmp/TIC1201247611_lc.png')

# >> pmin: minimum period (minutes)
# >> pmax: maximum period (days)
# >> qmin: minimum transit duraction (t_trans / period)
pmin=3
# pmax=True
pmax=0.05 # = 72 minutes
# qmin = 0.005
# qmax = 0.01
qmin=0.005
qmax=0.05
t, y, dy, period, bls_power_best, freqs, power = BLS(t,y,dy,pmin=pmin,pmax=pmax,
                                       qmin=qmin,qmax=qmax,remove=True)


plt.figure()
plt.plot(freqs, power, '.k', ms=1)
plt.xlim([0, 50])
plt.xlabel('Frequency')
plt.ylabel('BLS Power')
plt.savefig('/scratch/echickle/tmp/TIC1201247611_bls.png')

# >> compute phase curve
ts = TimeSeries(time=Time(t, format='jd'), data={'flux': y})
ts_folded = ts.fold(period=period*u.d) 

sorted_inds = np.argsort(ts_folded['time'])
folded_t = ts_folded['time'][sorted_inds]
folded_y = ts_folded['flux'][sorted_inds]

plt.figure()
plt.plot(folded_t.value, folded_y.value, '.k', ms=1)
plt.xlabel('Phase')
plt.ylabel('PDCSAP_FLUX')
plt.title('TIC 1201247611, pmin: '+str(qmin)+', pmax: '+str(pmax)+', qmin: '+str(qmin)+', qmax: '+str(qmax))
fname = '/scratch/echickle/tmp/TIC1201247611_pc_qmin'+str(qmin)+'_qmax'+str(qmax)+'_period'+str(period)+'_sig'+str(bls_power_best)+'.png'
plt.savefig(fname)
print('Saved '+fname)

# with open('/home/submit/echickle/foo.txt', 'w') as f:
#     f.write('Period '+str(period)+'\nPower '+str(bls_power_best))
# >> load white dwarf catalog
# wd_cat  = pd.read_csv('/home/echickle/work/tess_dwd/WDs.txt', header=None, sep='\s+')

# >> WDs observed in S30 short cadence
# s30 = np.loadtxt('/home/echickle/work/tess_dwd/all_targets_S030_v1.txt')
# n = len(np.intersect1d(s30[:,0], wd_cat[0]))

# >> load flux and time for white dwarfs in S30 (from TESS_Photometry.py)
# flux = np.load('/scratch/echickle/dwd/lc.npy').T
# time = np.load('/scratch/echickle/dwd/ts.npy') 
# ticid = np.load('/scratch/echickle/dwd/id.npy').astype('int')

# >> compute BLS
# for i in range(len(flux)):
#     y = flux[i]
#     inds = np.nonzero(~np.isnan(y))
#     t, y = time[inds], y[inds]
#     dy = np.ones(y.shape)
#     # default: pmin=3,pmax=True,qmin=2e-2,qmax=0.12
#     t, y, dy, period, bls_power_best = BLS(t,y,dy,pmin=3,pmax=True,qmin=2e-2,qmax=0.3,remove=True)

#     # >> compute phase curve
#     prefix = str(bls_power_best)+'-'+str(round(period,2))+'-TIC%016d'%ticid[i]+'_'
#     lc.make_phase_curve(t, y, period, output_dir=img_dir, prefix=prefix)

# pdb.set_trace()

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
