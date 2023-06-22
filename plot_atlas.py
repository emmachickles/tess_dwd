# ------------------------------------------------------------------------------

lc_dir = "/scratch/echickle/grnd_lc/"
atlas_dir = "/data/ATLAS/"
out_dir = "/scratch/echickle/vet_ATLAS/"
tess_dir = "/scratch/data/tess/lcur/ffi/"
gaia_tab = "/scratch/echickle/100pc_clean.fits"

bls=True

# ------------------------------------------------------------------------------

import lc_utils as lcu
import numpy as np
import sys
import os
import pdb
os.makedirs(out_dir, exist_ok=True)

# ------------------------------------------------------------------------------

# load data
gid, per, ra, dec = 282679289838317184

lcu.make_panel_plot(fname_atlas,fnames_ztf,tess_dir,ticid,cam,ccd,per,ra,dec,
                    per_atlas=per,
                    gaia_tab,out_dir,suffix,bls=bls)
