# output_dir = '/data/submit/echickle/out/'
# img_dir = output_dir + 's41-ztf/'
# data_dir = '/data/submit/echickle/data/sector_41_lc/'
output_dir = '/data/submit/tess/echickle/KBUCB-png/'
# img_dir = output_dir + 's0056-targ/'
data_dir = '/data/submit/tess/echickle/KBUCB-lc/'


# DWD: ./sector-05/tess2018319095959-s0005-0000000471013547-0125-s_lc.fits
# https://arxiv.org/pdf/2106.15104.pdf

import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt

import pdb
import os
import sys
import gc
import fnmatch
sys.path.insert(0, "/home/submit/echickle/work/")

from KBB_Utils.KBB_Utils.Period_Finding import BLS

from astropy.io import fits
from astropy.timeseries import TimeSeries
from astropy.time import Time
from astropy.coordinates import SkyCoord
import astropy.units as u
from wotan import flatten
import lc_utils as lcu
import gc

# >> ZTF
# fnames = os.listdir(data_dir)
# fnames = ['228.38793511208_70.62277160784.lc']
# fnames = ['266.405725_60.3234176.lc']
# ra = fnames[i].split('_')[0]
# dec = fnames[i].split('_')[1][:-3]

# >> AM CVn Gaia14aae
ticid = 1978047052

camccd = '-3-1'
while type(camccd) == type(None):
    for cam in [1,2,3,4]:
        for ccd in [1,2,3,4]:
            suffix = '-'+str(cam)+'-'+str(ccd)
            ticid_list = np.load(data_dir+'id'+suffix+'.npy').astype('int')
            ind = np.nonzero(ticid_list == ticid)
            if len(ind[0]) > 0:
                camccd = suffix

suffix = camccd
print(suffix)

ticid_list = np.load(data_dir+'id'+suffix+'.npy').astype('int')
ind = np.nonzero(ticid_list == ticid)[0][0]
print(ind)
t = np.load('/data/submit/echickle/data/sector_41_lc/ts'+suffix+'.npy')
y = np.load('/data/submit/echickle/data/sector_41_lc/lc'+suffix+'.npy')[ind]


for i in [5]:
    for j in [0.1]:

        # -- prep light curve (normalize, sigma-clip, detrend) -----------------

        t, y = lcu.prep_lc(t, y, n_std=i, wind=j)
        dy = np.ones(y.shape)

        # -- compute BLS -------------------------------------------------------
        #qmax=0.2
        t, y, dy, period, bls_power_best, freqs, power = \
            BLS(t,y,dy,pmin=7,pmax=0.25,qmin=0.005,qmax=0.2,remove=True)
            # BLS(t,y,dy,pmin=20,pmax=5,qmin=0.005,qmax=0.2,remove=True)
        
        print(period*1440)
        print(bls_power_best)

        # -- plot phase curve --------------------------------------------------
        # prefix = str(bls_power_best)+'-'+str(period*1440)+\
        #     '_ra'+ra+'_dec'+dec+'_'
        prefix = str(bls_power_best)+'-'+str(period*1440)+\
            '_TIC%016d'%ticid+'_'+camccd
            # '_ra'+ra+'_dec'+dec+'_'
        prefix = 'nstd'+str(i)+'-wind'+str(j)+'-' + prefix

        # inds = np.nonzero((freqs < 1440/20) * (freqs > 1440/30))
        # freqs, power = freqs[inds], power[inds]
        
        lcu.make_phase_curve(t, y, 24.76/1440, dy=dy, output_dir=img_dir,
                             prefix=prefix, freqs=freqs, power=power,
                             bins=100, ticid=ticid)

