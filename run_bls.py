# -- inputs --------------------------------------------------------------------

n_std = 5
detrend = "wotan"
wind = 0.1
pmin = 410 / 60 
pmax = 0.13 
qmin = 0.01
qmax = 0.15

sector, cam, ccd = 61, 3, 3

data_dir = "/scratch/data/tess/lcur/ffi/s%04d-lc/"%sector
output_dir = "/scratch/echickle/s%04d/"%sector
bls_dir    = output_dir + 's00{}-bls-{}-{}-230403/'.format(sector,cam,ccd)
# ls_dir     = output_dir + 's00{}-ls-{}-{}-230403/'.format(sector, cam,ccd)

# ------------------------------------------------------------------------------

import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt

import pdb
import os
import gc

from Period_Finding import BLS, LS_Astropy
import lc_utils as lcu

import sys

# ------------------------------------------------------------------------------

# cam = sys.argv[1]
# ccd = sys.argv[2]

os.makedirs(output_dir, exist_ok=True)
os.makedirs(bls_dir, exist_ok=True)
# os.makedirs(ls_dir, exist_ok=True)

fail_txt = output_dir + 'cam'+str(cam)+'-ccd'+str(ccd)+'-failed.txt'
with open(fail_txt, 'w') as f:
    f.write('')

# ------------------------------------------------------------------------------

suffix = '-'+str(cam)+'-'+str(ccd)+'.npy'
flux = np.load(data_dir+'lc'+suffix)
time = np.load(data_dir+'ts'+suffix)     
coord = np.load(data_dir+'co'+suffix)
ticid = np.load(data_dir+'id'+suffix).astype('int')

inds = np.argsort(time)
time, flux = time[inds], flux[:,inds]

# !! 
# ind = np.nonzero(ticid == 101433897) 
# flux = [flux[ind][0]]
# coord = [coord[ind][0]]
# ticid = [ticid[ind][0]]

# >> remove completed
fnames_ccd = os.listdir(bls_dir)
ticid_ccd = [int(f.split('_')[12][3:]) for f in fnames_ccd if f.split('.')[-1] == 'png']
ticid_ccd = np.array(ticid_ccd)
inter, comm1, comm2 = np.intersect1d(ticid, ticid_ccd, return_indices=True)
coord = np.delete(coord, comm1, axis=0)
flux = np.delete(flux, comm1, axis=0)
ticid = np.delete(ticid, comm1) 


# >> compute BLS
for i in range(len(flux)):
    if i % 50 == 0:
        print(str(i)+" / "+str(len(flux)))
    y = flux[i]
    t = time

    # -- prep light curve --------------------------------------------------
    t, y, flag = lcu.prep_lc(t, y, n_std=n_std, detrend=detrend, wind=wind)
    dy = np.ones(y.shape)
    
    if flag:
        with open(fail_txt, 'a') as f:
            f.write(str(ticid[i])+'\n')

    else:

        # -- compute BLS ---------------------------------------------------
        t, y, dy, period, bls_power_best, freqs, power, q, phi0 = \
            BLS(t,y,dy,pmin=pmin,pmax=pmax,qmin=qmin,qmax=qmax,remove=True)

        # -- plot phase curve ----------------------------------------------
        suffix = '_TIC%016d'%ticid[i]+'_s%04d_'%sector+'cam_'+\
                 str(cam)+'_ccd_'+str(ccd)+\
                 '_ra_{}_dec_{}_'.format(coord[i][0], coord[i][1])

        # !! save_npy
        lcu.vet_plot(t, y, freqs, power, q, phi0, output_dir=bls_dir,
                     ticid=ticid[i], suffix=suffix, bins=100, save_npy=False)

        # -- compute LS ----------------------------------------------------
        # _, _, _, ls_period, ls_power_best, ls_freqs, ls_power = \
        #     LS_Astropy(t,y,dy,pmax=pmax)
        
        # prefix = 'pow_'+str(ls_power_best)+'_per_'+str(round(ls_period*1440,5))+\
        #     '_TIC%016d'%ticid[i]+'_s%04d_'%sector+'cam_'+str(cam)+'_ccd_'+str(ccd)+\
        #     '_ra_{}_dec_{}_'.format(coord[i][0], coord[i][1])                

        # lcu.make_phase_curve(t, y, ls_period, output_dir=ls_dir,
        #                      prefix=prefix, freqs=ls_freqs, power=ls_power,
        #                      ticid=ticid[i], bins=100, bls=False, save_npy=True)

        
