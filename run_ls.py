# -- inputs --------------------------------------------------------------------

detrend = "wotan"
wind = 0.1
pmin = 410 / 60 
pmax = 0.13 
qmin = 0.01
qmax = 0.15

sector, cam, ccd = 57, 1, 4

import sys
# sector = int(sys.argv[1])
# cam    = int(sys.argv[2])
# ccd    = int(sys.argv[3])

wd_tab = "/home/echickle/data/WDs.txt"
wd_main = "/data/GaiaEDR3_WD_main.fits"
rp_ext = "/data/GaiaEDR3_WD_RPM_ext.fits"

data_dir = "/home/echickle/data/s%04d/"%sector+"/s%04d-lc/"%sector
output_dir = "/home/echickle/out/LS-s00{}-{}-{}/".format(sector, cam, ccd)

# ------------------------------------------------------------------------------

import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt

import pdb
import os
os.makedirs(output_dir, exist_ok=True)
import gc

from Period_Finding import LS_Astropy
import lc_utils as lcu

import sys

# ------------------------------------------------------------------------------

suffix = '-'+str(cam)+'-'+str(ccd)+'.npy'
flux = np.load(data_dir+'lc'+suffix)
time = np.load(data_dir+'ts'+suffix)     
coord = np.load(data_dir+'co'+suffix)
ticid = np.load(data_dir+'id'+suffix).astype('int')

inds = np.argsort(time)
time, flux = time[inds], flux[:,inds]

# !! 
ind = np.nonzero(ticid == 471015192) 
flux = [flux[ind][0]]
coord = [coord[ind][0]]
ticid = [ticid[ind][0]]

# # >> remove completed
# fnames_ccd = os.listdir(ls_dir)
# ticid_ccd = [int(f.split('_')[6][3:]) for f in fnames_ccd if f.split('.')[-1] == 'png']
# ticid_ccd = np.array(ticid_ccd)
# inter, comm1, comm2 = np.intersect1d(ticid, ticid_ccd, return_indices=True)
# coord = np.delete(coord, comm1, axis=0)
# flux = np.delete(flux, comm1, axis=0)
# ticid = np.delete(ticid, comm1) 

for i in range(len(flux)):
    if i % 50 == 0:
        print(str(i)+" / "+str(len(flux)))
    y = flux[i]
    t = time
    ra = coord[i][0]
    dec = coord[i][1]

    # -- prep light curve --------------------------------------------------
    t, y, flag = lcu.prep_lc(t, y, detrend=detrend, wind=wind)
    dy = np.ones(y.shape)*0.1

    # -- compute LS ----------------------------------------------------

    _, _, _, ls_period, ls_power_best, ls_freqs, ls_power = \
        LS_Astropy(t,y,dy,pmax=pmax)

    suffix='_TIC%016d'%ticid[i]+'_s%04d_'%sector+'cam_'+str(cam)+'_ccd_'\
        +str(ccd)+'_ra_{}_dec_{}'.format(ra, dec)    

    res=lcu.vet_plot(t, y, ls_freqs, ls_power,output_dir=output_dir,
                 objid=ticid[i], objid_type=None, suffix=suffix, ra=ra,
                 dec=dec, wd_main=wd_main, rp_ext=rp_ext, wd_tab=wd_tab,
                 bls=False)
    print(res)

    periods = [336.28, 336.68,  337.09, 379.55,  380.10, 380.66, 400.41, 400.99, 401.56, 448.79,  449.48, 450.16,  467.87, 468.69, 469.46,  562.70, 564.28, 565.59, 609.64, 611.15, 612.38]
    periods = np.array(periods)
    freqs = 86400 / periods

    fig, ax = plt.subplots(figsize=(10,4))
    for i in range(len(freqs)):
        ax.axvline(freqs[i], ls='dashed', color='r', alpha=0.3, lw=0.5)
    ax.axvline(1/ls_period, ls='dashed', color='b', alpha=0.2)
    ax.plot(ls_freqs, ls_power, '.k', ms=1)
    ax.set_xlabel('Frequency [1/days]')
    ax.set_ylabel('LS Power')
    fig.savefig(output_dir+'LS_periodogram'+suffix+'.png', dpi=300)
    print('Saved '+output_dir+'LS_periodogram'+suffix+'.png')

