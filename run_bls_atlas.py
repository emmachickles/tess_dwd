# -- inputs --------------------------------------------------------------------

pmin = 410 / 60 
pmax = 0.13 
qmin = 0.01
qmax = 0.15

data_dir = "/nobackup1c/users/echickle/ATLAS/"
output_dir = "/nobackup1c/users/echickle/out/"
bls_dir    = output_dir + 'bls'
ls_dir     = output_dir + 'ls'

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
os.makedirs(bls_dir, exist_ok=True)
os.makedirs(ls_dir, exist_ok=True)

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
for i in range(len(fnames)):
    if i % 50 == 0:
        print(str(i)+" / "+str(len(flux)))
    y = flux[i]
    t = time

    t, y, dy = load_atlas_lc(fname_atlas)
    t, y, dy, period, bls_power_best, freqs, power, q, phi0 = \
        BLS(t,y,dy,pmin=pmin,pmax=pmax,qmin=qmin,qmax=qmax,remove=True)
    prefix1 = 'ATLAS_'+prefix+'_'
    vet_plot(t, y, freqs, power, q, phi0, output_dir=out_dir+prefix1,
             ticid=ticid, dy=dy)

    

    # -- compute LS ----------------------------------------------------

    _, _, _, ls_period, ls_power_best, ls_freqs, ls_power = \
        LS_Astropy(t,y,dy,pmax=pmax)

    suffix='_TIC%016d'%ticid[i]+'_s%04d_'%sector+'cam_'+str(cam)+'_ccd_'\
        +str(ccd)+'_ra_{}_dec_{}_'.format(coord[i][0], coord[i][1])    

    lcu.vet_plot(t, y, ls_freqs, ls_power,output_dir=ls_dir,
                 ticid=ticid[i], suffix=suffix, bins=100, save_npy=False,
                 bls=False)

