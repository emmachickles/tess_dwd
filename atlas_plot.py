# ------------------------------------------------------------------------------

lc_dir = "/scratch/echickle/grnd_lc/"

lc_f = lc_dir + "LDSS_230421.txt" # !!
sector, cam = 62, 3 # !!

out_dir = "/scratch/echickle/LDSS_230421_ATLAS/"
# out_dir = "/scratch/echickle/s%04d/"%sector \
#           +"s%04d-"%sector+str(cam)+"-crossmatch/"
tess_dir = "/scratch/data/tess/lcur/ffi/s%04d-lc/"%sector

gaia_tab = "/scratch/echickle/100pc_clean.fits"
bls=True

# ------------------------------------------------------------------------------

import lc_utils as lcu
import numpy as np
import sys
import os
import pdb

# ------------------------------------------------------------------------------

lc_tab = np.loadtxt(lc_f, delimiter=",", dtype="str", skiprows=1)
if len(lc_tab.shape) == 1:
    lc_tab = np.expand_dims(lc_tab, 0)
ticid = np.unique(lc_tab[:,0])
os.makedirs(out_dir, exist_ok=True)

pmin = 410 / 60 
pmax = 0.13 
qmin = 0.01
qmax = 0.15

# >> remove completed
# fnames = os.listdir(out_dir)
# ticid_out = [int(f.split('_')[13][3:]) for f in fnames if f.split('_')[0] == 'TESS']
# ticid_out = np.array(ticid_out)
# inter, comm1, comm2 = np.intersect1d(ticid, ticid_out, return_indices=True)
# ticid = np.delete(ticid, comm1) 

# !!
ticid = ['808364853']

for i in range(len(ticid)):
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
                        out_dir,bls=bls,pmin=pmin,pmax=pmax,qmin=qmin,qmax=qmax)
