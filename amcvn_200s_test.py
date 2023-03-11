# Testing q factor range
# AM CVn TIC 1201247611

import os
import pdb
import numpy as np
import lc_utils as lcu
from Period_Finding import BLS, LS_Astropy

ticid = 1201247611
data_dir = "/scratch/data/tess/lcur/ffi/s0057-lc/"

# Where to save BLS, LS plots 
output_dir = "/scratch/echickle/s0057/AMCVn_test/"
bls_dir    = output_dir+"bls/"
ls_dir     = output_dir+"ls/"
os.makedirs(bls_dir, exist_ok=True)
os.makedirs(ls_dir, exist_ok=True)

# Light curve extraction parameters
n_std = 5
detrend = "polyfit"

# BLS & LS period minimum & maximum
pmin = 410 / 60
pmax = 0.13
qmin = 0.01
qmax = 0.15

# ------------------------------------------------------------------------------

# # Make text file 
# ft = open(output_dir+"out_table.txt", "w")
# cols = ["wind", "per[min]", "sig"]
# ft.write(("{}\t\t"*len(cols)).format(*cols)+"\n")
# ft.close()

# Get time and flux array
ticid_list = np.load(data_dir+'id-4-3.npy')
ind = np.nonzero(ticid_list == ticid)[0][0]
t = np.load(data_dir+'ts-4-3.npy')
y = np.load(data_dir+'lc-4-3.npy')[ind]
coord = np.load(data_dir+'co-4-3.npy')[ind]

# qrng = [[0.005, 0.05], [0.01, 0.05], [0.05,0.3]]
# wrng = ["poly8"]


t, y, flag = lcu.prep_lc(t, y, n_std=n_std, detrend=detrend)
dy = np.ones(y.shape) * 0.1

# Run BLS
_, _, _, per, pow_best, bls_frq, bls_pow, dur, epo, delta = \
        BLS(t,y,dy,pmin=pmin,pmax=pmax,qmin=qmin,qmax=qmax,remove=True)
_, _, _, ls_per, ls_pow_best, ls_frq, ls_pow = LS_Astropy(t,y,dy,pmax=pmax)

prefix = 'pow_'+str(pow_best)+'_delta_'+str(round(delta,5))+\
         '_per_'+str(round(per*1440,5))+\
    '_TIC%016d'%ticid+'_dur_'+str(dur)+'_epo_'+str(epo)+\
    '_ra_{}_dec_{}_'.format(coord[0], coord[1])                

lcu.make_phase_curve(t, y, per, output_dir=bls_dir, prefix=prefix,
                     freqs=bls_frq, power=bls_pow, ticid=ticid, bins=100,
                     dur=dur, epo=epo)

prefix = 'pow_'+str(ls_pow_best)+'_per_'+str(round(ls_per*1440,5))+\
    '_TIC%016d'%ticid+'_ra_{}_dec_{}_'.format(coord[0], coord[1])                

lcu.make_phase_curve(t, y, ls_per, output_dir=ls_dir, prefix=prefix,
                     freqs=ls_frq, power=ls_pow, ticid=ticid, bins=100,
                     bls=False)


# ft = open(output_dir+"out_table.txt", "a")
# ft.write(("{}\t\t"*len(cols)).format(wind, round(per*1440,5),
#                                      round(pow_best,5))+"\n")
# ft.close()

