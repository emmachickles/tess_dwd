# -- inputs --------------------------------------------------------------------

n_std = 5
detrend = "wotan"
wind = 0.1
pmin = 400 / 60 
pmax = 10
qmin = 0.01
qmax = 0.15

# output_dir = "/scratch/echickle/BLS_Test_Nonastro/"
# ticid  = [836153158, 830748654, 778219027, 154178878, 803454388, 803568171]
# sector = [61,        61,        61,        61,        61,        61]
# cam    = [2,         2,         2,         2,         1,         1]
# ccd    = [3,         3,         1,         2,         1,         2]

# output_dir = "/scratch/echickle/BLS_Test_Real/"
# ticid  = [803489769, 36085812, 452954413, 455206965, 800042858, 1201247611]
# sector = [61,        61,       61,        61,        61,        57]
# cam    = [1,         1,        1,         1,         1,         4]
# ccd    = [1,         1,        4,         4,         3,         3]

output_dir = "/scratch/echickle/BLS_Phase_Entropy/"
ticid  = [803280089, 803332316, 802868096, 803849868, 803369981, 803313532]
sector = [61,        61,        61,        61,        61,        61]
cam    = [1,         1,         1,         1,         1,         1]
ccd    = [1,         1,         1,         1,         1,         1]

# output_dir = "/scratch/echickle/MGAB/"
# ticid  = [803489769, 36085812,  835734923, 875850017, 471014834]
# sector = [61,        61,        61,        62,        62]
# cam    = [1,         1,         2,         1,         1]
# ccd    = [1,         1,         4,         2,         3]

# output_dir = "/scratch/echickle/KB_UCB/" # incomplete! 
# ticid  = [2040677137, 746047312]
# sector = [57,         61]
# cam    = [2,          2]
# ccd    = [3,          4]

# output_dir = "/scratch/echickle/LDSS_230421/"
# ticid  = [767706310, 875850017, 826164830, 193092806, 808364853]
# sector = [61,        62,        62,        62,        62]
# cam    = [3,         1,         2,         2,         3]
# ccd    = [3,         2,         1,         3,         4]

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

os.makedirs(output_dir, exist_ok=True)
bls_dir    = output_dir + 'bls/'
ls_dir     = output_dir + 'ls/'
os.makedirs(bls_dir, exist_ok=True)
os.makedirs(ls_dir, exist_ok=True)

# ------------------------------------------------------------------------------

result = []
# >> compute BLS
for i in range(len(ticid)):
    print(ticid[i])

    data_dir = "/scratch/data/tess/lcur/ffi/s%04d-lc/"%sector[i]

    suffix = '-'+str(cam[i])+'-'+str(ccd[i])+'.npy'

    ticid_list = np.load(data_dir+'id'+suffix).astype('int')
    ind = np.nonzero(ticid_list == ticid[i])[0][0]
    y = np.load(data_dir+'lc'+suffix)[ind]
    coord = np.load(data_dir+'co'+suffix)[ind]

    t = np.load(data_dir+'ts'+suffix)     

    # # !! 
    # from astropy.io import fits
    # hdul = fits.open("/home/echickle/hlsp_tglc_tess_ffi_gaiaid-1916879054217244544-s0057-cam2-ccd3_tess_v1_llc.fits")
    # t = hdul[1].data['time']
    # y = hdul[1].data['cal_psf_flux']


    # # # !!
    # inds = np.nonzero( t> 59870 )
    # t, y = t[inds], y[inds]

    # -- prep light curve --------------------------------------------------
    # if i == 1: pdb.set_trace()
    t, y, flag = lcu.prep_lc(t, y, n_std=n_std, detrend=detrend, wind=wind)
    dy = np.ones(y.shape)*0.1
    
    # -- compute BLS ---------------------------------------------------
    t, y, dy, period, bls_power_best, freqs, power, q, phi0 = \
        BLS(t,y,dy,pmin=pmin,pmax=pmax,qmin=qmin,qmax=qmax,remove=True)

    # -- plot phase curve ----------------------------------------------
    suffix = '_TIC%016d'%ticid[i]+'_s%04d_'%sector[i]+'cam_'+\
             str(cam[i])+'_ccd_'+str(ccd[i])+\
             '_ra_{}_dec_{}_'.format(coord[0], coord[1])

    res = lcu.vet_plot(t, y, freqs, power, q, phi0, output_dir=bls_dir,
                 objid=ticid[i], suffix=suffix, bins=100, save_npy=False)
    res = [ticid[i], coord[0], coord[1]] + list(res)
    result.append(res)

# ticid, ra, dec, sig, snr, wid, period, period_min, q, phi0, dur, epo, rp, nt, dphi
np.savetxt(output_dir+'GPU.result', np.array(result),
           fmt='%s,%10.5f,%10.5f,%10.5f,%10.5f,%i,%10.5f,%10.5f,%10.5f,%10.5f,%10.5f,%10.5f,%10.5f,%i,%10.5f',
           header='ticid, ra, dec, sig, snr, wid, period, period_min, q, phi0, dur, epo, rp, nt, dphi')
print('Saved '+output_dir+'GPU.result')

    # # period = 0.038365738 # TIC 2040677137
    # # period = 0.038365738/2 # TIC 2040677137
    # peak = np.argmin(np.abs(freqs - 1/period) )
    # # # peak = np.argmax(power)
    # nearpeak = 3000
    # fig, ax = plt.subplots()
    # ax.axvline(1/period, color='r')
    # ax.plot(freqs[max(0,peak-nearpeak):peak+nearpeak],
    #            power[max(0,peak-nearpeak):peak+nearpeak], '.k', ms=1)
    # ax.set_xlabel('Frequency [1/days]') 
    # ax.set_ylabel('BLS Power')
    # fig.savefig(output_dir+'TIC%016d'%ticid[i]+'_zoom_bls.png')
    # print(output_dir+'TIC%016d'%ticid[i]+'_zoom_bls.png')
    # # pdb.set_trace()

    # # -- compute LS ----------------------------------------------------
    # _, _, _, ls_period, ls_power_best, ls_freqs, ls_power = \
    #     LS_Astropy(t,y,dy,pmax=pmax)

    # suffix='_TIC%016d'%ticid[i]+'_s%04d_'%sector[i]+'cam_'+str(cam[i])+'_ccd_'+str(ccd[i])+\
    #     '_ra_{}_dec_{}_'.format(coord[0], coord[1])                

    # lcu.vet_plot(t, y, ls_freqs, ls_power,output_dir=ls_dir,
    #              ticid=ticid[i], suffix=suffix, bins=100, save_npy=False,
    #              bls=False)

    # period = 0.038365738/2 # TIC 2040677137                                      
    # peak = np.argmin(np.abs(ls_freqs - 1/period) )
    # # # peak = np.argmax(power)                                                  
    # nearpeak = 3000
    # fig, ax = plt.subplots()
    # ax.axvline(1/period, color='r')
    # ax.plot(ls_freqs[max(0,peak-nearpeak):peak+nearpeak],
    #         ls_power[max(0,peak-nearpeak):peak+nearpeak], '.k', ms=1)
    # ax.set_xlabel('Frequency [1/days]')
    # ax.set_ylabel('LS Power')
    # fig.savefig(output_dir+'TIC%016d'%ticid[i]+'_zoom_ls.png')
    # print(output_dir+'TIC%016d'%ticid[i]+'_zoom_ls.png')

    # pdb.set_trace()
