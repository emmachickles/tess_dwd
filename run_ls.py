import os
import sys

# sector = int(sys.argv[1])
# cam    = int(sys.argv[2])
# ccd    = int(sys.argv[3])

# -- inputs --------------------------------------------------------------------

detrend = "wotan"
wind = 0.1
pmin = 400 / 60 
pmax = 10
n_std=5

wid_threshold=6
pow_threshold=25

# >>ENGAGING
wd_tab= "/nobackup1c/users/echickle/WDs.txt"
wd_main = "/nobackup1c/users/echickle/GaiaEDR3_WD_main.fits"
rp_ext = "/nobackup1c/users/echickle/GaiaEDR3_WD_RPM_ext.fits"
qflag_dir = "/nobackup1c/users/echickle/QLPqflags/"

# >> UZAY

# sector, cam, ccd = 65, 1, 2

# wd_tab = "/scratch/echickle/WDs.txt"
# wd_main = "/scratch/echickle/GaiaEDR3_WD_main.fits"
# rp_ext = "/scratch/echickle/GaiaEDR3_WD_RPM_ext.fits"
# qflag_dir = "/scratch/echickle/QLPqflags/"
# data_dir = "/scratch/data/tess/lcur/ffi/s%04d-lc/"%sector

# import os
# os.makedirs("/scratch/echickle/s%04d_LS/"%sector, exist_ok=True)
# output_dir = "/scratch/echickle/s%04d_LS/"%sector+"cam{}-ccd{}/".format(cam,ccd)
# stat_dir = "/scratch/echickle/s%04d_LS/"%sector+"stats/"

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

# os.makedirs(output_dir, exist_ok=True)
# os.makedirs(stat_dir, exist_ok=True)

# ------------------------------------------------------------------------------

# suffix = '-'+str(cam)+'-'+str(ccd)+'.npy'
# flux = np.load(data_dir+'lc'+suffix)
# time = np.load(data_dir+'ts'+suffix)     
# coord = np.load(data_dir+'co'+suffix)
# cn = np.load(data_dir+'cn'+suffix)
# ticid = np.load(data_dir+'id'+suffix).astype('int')

# inds = np.argsort(time)
# time, flux = time[inds], flux[:,inds]

# !! 
# ind = np.nonzero(ticid == 150808542) 
# flux = [flux[ind][0]]
# coord = [coord[ind][0]]
# ticid = [ticid[ind][0]]

# # >> remove completed
# fnames_ccd = os.listdir(ls_dir)
# ticid_ccd = [int(f.split('_')[6][3:]) for f in fnames_ccd if f.split('.')[-1] == 'png']
# ticid_ccd = np.array(ticid_ccd)
# inter, comm1, comm2 = np.intersect1d(ticid, ticid_ccd, return_indices=True)
# coord = np.delete(coord, comm1, axis=0)
# flux = np.delete(flux, comm1, axis=0)
# ticid = np.delete(ticid, comm1) 

# result_list = []

def run_process(p):
    sector, cam, ccd, ticid, data_dir, ls_dir = p

    print('Starting S{}-{}-{} TIC{}'.format(sector, cam, ccd, ticid))

    # >> load data
    ticid = np.int64(ticid)
    suffix = "-{}-{}.npy".format(cam, ccd)
    y = np.load(data_dir+'lc'+suffix)
    cn = np.load(data_dir+'cn'+suffix)
    # y = np.load(data_dir+'lc-2-2-ap1.1-in1.8-out2.3.npy')
    t = np.load(data_dir+'ts'+suffix)
    coord = np.load(data_dir+'co'+suffix)
    ticid_ccd = np.load(data_dir+'id'+suffix).astype('int')
    ind = np.nonzero(ticid_ccd == ticid)[0][0]
    y, coord = y[ind], coord[ind]
    ra, dec = coord[0], coord[1]
    print('Loaded S{}-{}-{} TIC{}'.format(sector, cam, ccd, ticid))

    # >> detrend and remove outliers
    t, y, cn = lcu.rm_qflag(t, y, cn, qflag_dir, sector, cam, ccd)
    t, y, flag = lcu.prep_lc(t, y, n_std=n_std, detrend=detrend, wind=wind)

    dy = np.ones(y.shape)
    freqs_to_remove=lcu.rm_freq_tess()
    _, _, _, ls_period, ls_freqs, ls_power = LS(t,y,dy,freqs_to_remove=freqs_to_remove, pmin=200/60.)

    suffix='_TIC%016d'%ticid+'_s%04d_'%sector+'cam_'+str(cam)+'_ccd_'\
        +str(ccd)+'_ra_{}_dec_{}'.format(ra, dec)    

    res=lcu.vet_plot(t, y, ls_freqs, ls_power, output_dir=ls_dir,
                 objid=ticid, objid_type='TICID', suffix=suffix, ra=ra,
                 dec=dec, wd_main=wd_main, rp_ext=rp_ext, wd_tab=wd_tab,
                     bls=False, wid_threshold=wid_threshold, pow_threshold=pow_threshold)
    # print(res) # sig, wid, period, period_min, dphi
    return [np.int64(ticid), ra, dec] + list(res)

    # -- true frequencies ----------------------------------------
    # freqs = np.array([1440./55.24666272])


    # # PG 0122+200
    # periods = [609.5729, 467.8656, 450.1753, 449.5010, 401.5558, 400.9785, 400.4042, 380.1111, 336.2804]
    # periods = np.array(periods)
    # freqs = 86400 / periods

    # # GD 29-38
    # freqs = np.array([0.9847, 1.0742, 1.2126, 1.4771, 1.6317]) # hm check
    # freqs = freqs * 1e3

    # # WD J1527
    # periods = [745.87715, 704.24127, 702.89501, 701.75412, 649.39711, 351.46205]
    # periods = np.array(periods)
    # freqs = 86400 / periods

    # fig, ax = plt.subplots(figsize=(10,4))
    # for i in range(len(freqs)):
    #     ax.axvline(freqs[i], ls='dashed', color='r', alpha=0.3, lw=0.5)
    # ax.axvline(1/ls_period, ls='dashed', color='b', alpha=0.2)
    # ax.plot(ls_freqs, ls_power, '.k', ms=1)
    # ax.set_xlabel('Frequency [1/days]')
    # ax.set_ylabel('LS Power')
    # fig.savefig(output_dir+'LS_periodogram'+suffix+'.png', dpi=300)
    # print('Saved '+output_dir+'LS_periodogram'+suffix+'.png')
