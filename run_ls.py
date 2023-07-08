# -- inputs --------------------------------------------------------------------

detrend = "wotan"
wind = 0.1
pmin = 410 / 60 
pmax = 0.13 
qmin = 0.01
qmax = 0.15

sector, cam, ccd = 65, 1, 2

import sys
# sector = int(sys.argv[1])
# cam    = int(sys.argv[2])
# ccd    = int(sys.argv[3])

# wd_tab = "/home/echickle/data/WDs.txt"
# wd_main = "/data/GaiaEDR3_WD_main.fits"
# rp_ext = "/data/GaiaEDR3_WD_RPM_ext.fits"
# data_dir = "/home/echickle/data/s%04d/"%sector+"/s%04d-lc/"%sector
# output_dir = "/home/echickle/out/LS-s00{}-{}-{}/".format(sector, cam, ccd)

# >> ENGAGING
wd_tab= "/nobackup1c/users/echickle/WDs.txt"
wd_main = "/nobackup1c/users/echickle/GaiaEDR3_WD_main.fits"
rp_ext = "/nobackup1c/users/echickle/GaiaEDR3_WD_RPM_ext.fits"
qflag_dir = "/nobackup1c/users/echickle/QLPqflags/"
data_dir = "/nobackup1c/users/echickle/TESS_Lightcurves/s%04d-lc/"%sector

# >> UZAY
# wd_tab = "/scratch/echickle/WDs.txt"
# wd_main = "/scratch/echickle/GaiaEDR3_WD_main.fits"
# rp_ext = "/scratch/echickle/GaiaEDR3_WD_RPM_ext.fits"
# qflag_dir = "/scratch/echickle/QLPqflags/"
# data_dir = "/scratch/data/tess/lcur/ffi/s%04d-lc/"%sector


# ------------------------------------------------------------------------------

from Period_Finding import LS_Astropy, LS
import lc_utils as lcu
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import pdb
import gc
import sys

import os
os.makedirs("/scratch/echickle/s%04d_LS/"%sector, exist_ok=True)
output_dir = "/scratch/echickle/s%04d_LS/"%sector+"cam{}-ccd{}/".format(cam,ccd)
stat_dir = "/scratch/echickle/s%04d_LS/"%sector+"stats/"
os.makedirs(output_dir, exist_ok=True)
os.makedirs(stat_dir, exist_ok=True)

# ------------------------------------------------------------------------------

suffix = '-'+str(cam)+'-'+str(ccd)+'.npy'
flux = np.load(data_dir+'lc'+suffix)
time = np.load(data_dir+'ts'+suffix)     
coord = np.load(data_dir+'co'+suffix)
ticid = np.load(data_dir+'id'+suffix).astype('int')

inds = np.argsort(time)
time, flux = time[inds], flux[:,inds]

# !! 
ind = np.nonzero(ticid == 150808542) 
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

result_list = []

for i in range(len(flux)):
    if i % 50 == 0:
        print(str(i)+" / "+str(len(flux)))
    y = flux[i]
    t = time
    ra = coord[i][0]
    dec = coord[i][1]

    # -- prep light curve --------------------------------------------------
    t, y, cn = lcu.rm_qflag(t, y, cn, qflag_dir, sector, cam, ccd)
    t, y, flag = lcu.prep_lc(t, y, detrend=detrend, wind=wind)
    dy = np.ones(y.shape)*0.1

    # -- compute LS ----------------------------------------------------

    # _, _, _, ls_period, ls_power_best, ls_freqs, ls_power = \
    #     LS_Astropy(t,y,dy,pmax=pmax)
    freqs_to_remove=lcu.rm_freq_tess()
    _, _, _, ls_period, ls_freqs, ls_power = LS(t,y,dy,freqs_to_remove=freqs_to_remove, pmin=200/60.)

    suffix='_TIC%016d'%ticid[i]+'_s%04d_'%sector+'cam_'+str(cam)+'_ccd_'\
        +str(ccd)+'_ra_{}_dec_{}'.format(ra, dec)    

    res=lcu.vet_plot(t, y, ls_freqs, ls_power,output_dir=output_dir,
                 objid=ticid[i], objid_type=None, suffix=suffix, ra=ra,
                 dec=dec, wd_main=wd_main, rp_ext=rp_ext, wd_tab=wd_tab,
                 bls=False)
    print(res) # sig, wid, period, period_min, dphi
    result_list.append([np.int64(ticid[i]), ra, dec] + list(res))

    # -- true frequencies ----------------------------------------
    # freqs = np.array([1440./55.24666272])


    # PG 0122+200
    periods = [609.5729, 467.8656, 450.1753, 449.5010, 401.5558, 400.9785, 400.4042, 380.1111, 336.2804]
    periods = np.array(periods)
    freqs = 86400 / periods

    # GD 29-38
    freqs = np.array([0.9847, 1.0742, 1.2126, 1.4771, 1.6317]) # hm check
    freqs = freqs * 1e3

    # WD J1527
    periods = [745.87715, 704.24127, 702.89501, 701.75412, 649.39711, 351.46205]
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

np.savetxt(stat_dir+'stat-{}-{}-{}.result'.format(sector, cam, ccd), np.array(result_list),
           fmt='%s,%10.7f,%10.7f,%10.5f,%i,%10.5f,%10.8f,%10.5f',
           header='ticid, ra, dec, sig, wid, per, per_min, dphi')
