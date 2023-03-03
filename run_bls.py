# -- inputs --------------------------------------------------------------------

n_std = 3
wind  = 0.1
pmin = 400 / 60 # Nyquist in minutes
pmax = 0.15 # >> 3 hours = 0.125 days
qmin = 0.01
qmax = 0.05

# pmax = 0.25
# qmin = 0.005
# qmax = 0.05


# output_dir = '/home/echickle/out/'
# data_dir   = '/home/echickle/data/s0061/s0061-lc/'

# bls_dir    = output_dir + 's0061-bls-230228/'
# ls_dir    = output_dir + 's0061-ls-230228/'
# diag_dir   = output_dir + 's0061-diag-230228/'

import sys
cam = sys.argv[1]
ccd = sys.argv[2]

output_dir = '/data/submit/tess/echickle/'
data_dir   = '/data/submit/echickle/s0061-lc/'

bls_dir    = output_dir + 's0061-bls-{}-{}-230302/'.format(cam,ccd)
ls_dir     = output_dir + 's0061-ls-{}-{}-230302/'.format(cam,ccd)
diag_dir    = output_dir + 's0061-diag-{}-{}-230302/'.format(cam,ccd)

# ------------------------------------------------------------------------------

import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt

import pdb
import os
import gc
import fnmatch
sys.path.insert(0, "/home/submit/echickle/work/")

from Period_Finding import BLS, LS_Full
# from KBB_Utils.KBB_Utils.Period_Finding import BLS, LS_Full

from astropy.io import fits
from astropy.timeseries import TimeSeries
from astropy.time import Time
from astropy.coordinates import SkyCoord
import astropy.units as u
from wotan import flatten
import lc_utils as lcu
import gc


os.makedirs(bls_dir, exist_ok=True)
os.makedirs(ls_dir, exist_ok=True)
os.makedirs(diag_dir, exist_ok=True)

fail_txt = output_dir + 'cam'+str(cam)+'-ccd'+str(ccd)+'-failed.txt'
with open(fail_txt, 'w') as f:
    f.write('')

suffix = '-'+str(cam)+'-'+str(ccd)+'.npy'
flux = np.load(data_dir+'lc'+suffix)
time = np.load(data_dir+'ts'+suffix)     
coord = np.load(data_dir+'co'+suffix)
ticid = np.load(data_dir+'id'+suffix).astype('int')

inds = np.argsort(time)
time, flux = time[inds], flux[:,inds]

# ind = np.nonzero(ticid == 832987204)
# flux = [flux[ind][0]]
# coord = [coord[ind][0]]
# ticid = [ticid[ind][0]]

# >> remove completed
fnames_ccd = os.listdir(bls_dir)
ticid_ccd = [int(f.split('_')[4][3:]) for f in fnames_ccd]
ticid_ccd = np.array(ticid_ccd)
inter, comm1, comm2 = np.intersect1d(ticid, ticid_ccd, return_indices=True)
coord = np.delete(coord, comm1, axis=0)
flux = np.delete(flux, comm1, axis=0)
ticid = np.delete(ticid, comm1) 

# with open('/home/submit/echickle/out_cam{}_ccd{}.txt'.format(cam, ccd), 'w') as f:
#     f.write(str(len(flux))+'\n')

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

        # with open('/home/submit/echickle/out_cam{}_ccd{}.txt'.format(cam, ccd), 'a') as f:            
        #     f.write('i='+str(i)+' '+str(ticid[i])+' failed \n')            
    else:

        # -- compute BLS ---------------------------------------------------
        _, _, _, period, bls_power_best, freqs, power, dur, epo = \
            BLS(t,y,dy,pmin=pmin,pmax=pmax,qmin=qmin,qmax=qmax,remove=True)

        # pdb.set_trace()

        # # -- compute LS --------------------------------------------------
        # _, _, _, ls_period, ls_power_best, ls_freqs, ls_power = LS_Full(t,y,dy,pmin=pmin,
        #                                                        pmax=0.25)

        # -- plot phase curve ----------------------------------------------
        prefix = 'pow_'+str(bls_power_best)+'_per_'+str(round(period*1440,5))+\
            '_TIC%016d'%ticid[i]+'_cam_'+str(cam)+'_ccd_'+str(ccd)+\
            '_dur_'+str(dur)+'_epo_'+str(epo)+\
            '_ra_{}_dec_{}_'.format(coord[i][0], coord[i][1])                

        # lcu.plot_phase_curve(t, y, period, bls_dir, bins=100, prefix=prefix)
        lcu.make_phase_curve(t, y, period, dy=dy, output_dir=bls_dir,
                             prefix=prefix, freqs=freqs, power=power,
                             ticid=ticid[i], bins=100)

        # with open('/home/submit/echickle/out_cam{}_ccd{}.txt'.format(cam, ccd), 'a') as f:             
        #     f.write('i='+str(i)+' '+str(ticid[i])+' success \n')   
        
        # prefix = 'pow_'+str(ls_power_best)+'_per_'+str(round(ls_period*1440,5))+\
        #     '_TIC%016d'%ticid[i]+'_cam_'+str(cam)+'_ccd_'+str(ccd)+\
        #     '_ra_{}_dec_{}_'.format(coord[i][0], coord[i][1])                

        # lcu.make_phase_curve(t, y, ls_period, dy=dy, output_dir=ls_dir,
        #                      prefix=prefix, freqs=ls_freqs, power=ls_power,
        #                      ticid=ticid[i], bins=100, bls=False)
        
# with open('/home/submit/echickle/out_cam{}_ccd{}.txt'.format(cam, ccd), 'a') as f:      
#     f.write('Done!')   
        
# # >> diagnostic plot
# fnames = os.listdir(img_dir)
# period = []
# for i in range(len(fnames)):
#     period.append(float(fnames[i].split('_')[3]))
# period = np.array(period)
# plt.figure()
# plt.hist(period, bins=200)
# plt.xlabel('Period (minutes)')
# plt.ylabel('Num targets in S0061 Cam {} CCD {}'.format(cam, ccd))
# plt.yscale('log')
# plt.savefig(diag_dir+'hist_period_nstd_{}_wind_{}.png'.format(n_std, wind))

# plt.figure()
# plt.hist(1440/period, bins=200)
# plt.xlabel('Frequency (1/days)')
# plt.ylabel('Num targets in S0061 Cam {} CCD {}'.format(cam, ccd))
# plt.yscale('log')
# plt.savefig(diag_dir+'hist_freq_nstd_{}_wind_{}.png'.format(n_std, wind))



        # if len(fname.split('-')) > 3:
        #     N_ap = fname.split('-')[3]
        #     N_in = 'in'+fname.split('-')[4][4:7]
        #     N_out = 'out'+fname.split('-')[4][9:12]

        #     prefix = 'pow_'+str(bls_power_best)+'_per_'+str(round(period*1440,5))+\
        #         '_TIC%016d'%ticid[i]+'_cam_'+str(cam)+'_ccd_'+str(ccd)+\
        #         '_dur_'+str(dur)+'_epo_'+str(epo)+\
        #         '_Nap_'+N_ap+'_Nin_'+N_in+'_Not_'+N_out+\
        #         '_ra_{}_dec_{}_'.format(coord[i][0], coord[i][1])
        # else:
