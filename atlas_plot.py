# ------------------------------------------------------------------------------


lc_dir = "/scratch/echickle/ground_lc/"
out_dir = "/home/echickle/data/s0061/s0061-bls-1-3-230322-plot/" 
ls = False

# ------------------------------------------------------------------------------

import lc_utils as lcu
import sys
import os

# ------------------------------------------------------------------------------

gaia_tab = "/home/echickle/data/100pc_clean.fits"
wd_tab = "/home/echickle/data/WDs.txt"
tess_dir = "/home/echickle/data/s0061/s0061-lc/"
atlas_dir = "/data/ATLAS/"
os.makedirs(out_dir, exist_ok=True)

pmin = 410 / 60 
pmax = 0.13 
qmin = 0.01
qmax = 0.15

for fname in fnames:

    lcu.make_panel_plot(fname_tess,fname_atlas,fnames_ztf,tess_dir,gaia_tab,
                        out_dir,ls=ls)
