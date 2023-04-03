# ------------------------------------------------------------------------------


lc_dir = "/scratch/echickle/grnd_lc/"
# lc_f = lc_dir + "s0061-test.txt"
lc_f = lc_dir + "s0061-cam2-ccd3.txt"

# out_dir = "/home/echickle/data/s0061/s0061-test/" 
out_dir = "/home/echickle/data/s0061/s0061-cam2-ccd4/" 

gaia_tab = "/scratch/echickle/100pc_clean.fits"
tess_dir = "/scratch/data/tess/lcur/ffi/s0061-lc/"
ls = False

# ------------------------------------------------------------------------------

import lc_utils as lcu
import numpy as np
import sys
import os
import pdb

# ------------------------------------------------------------------------------

lc_tab = np.loadtxt(lc_f, delimiter=",", dtype="str", skiprows=1)
ticid = np.unique(lc_tab[:,0])

os.makedirs(out_dir, exist_ok=True)

pmin = 410 / 60 
pmax = 0.13 
qmin = 0.01
qmax = 0.15

for i in range(0, len(ticid)):
    print(str(i)+'/'+str(len(ticid)))
    lc_info = lc_tab[lc_tab[:,0] == ticid[i]]
    fname_tess = lc_info[lc_info[:,1] == 'TESS'][0][2]
    fname_atlas = lc_info[lc_info[:,1] == 'ATLAS'][:,2]
    if len(fname_atlas) > 0:
        fname_atlas = lc_dir+fname_atlas[0]
    else:
        fname_atlas = None
    fnames_ztf = list(lc_info[lc_info[:,1] == 'ZTF'][:,2])
    for j in range(len(fnames_ztf)):
        fnames_ztf[j] = lc_dir+fnames_ztf[j]

    lcu.make_panel_plot(fname_tess,fname_atlas,fnames_ztf,tess_dir,gaia_tab,
                        out_dir,ls=ls,pmin=pmin,pmax=pmax,qmin=qmin,qmax=qmax)
