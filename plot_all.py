# ------------------------------------------------------------------------------

atlas_dir = "/data/ATLAS/"
tess_dir = "/home/echickle/data/"
out_dir = "/home/echickle/out/vet/"
gaia_tab = "/home/echickle/data/100pc_clean.fits"
wd_tab = "/home/echickle/data/WDs.txt"
wd_main = "/data/GaiaEDR3_WD_main.fits"
rp_ext = "/data/GaiaEDR3_WD_RPM_ext.fits"

bls=True
clip=False

# ------------------------------------------------------------------------------

import lc_utils as lcu
import numpy as np
import sys
import os
import pdb
os.makedirs(out_dir, exist_ok=True)

# ------------------------------------------------------------------------------

# load data
ra_list = [294.59586]
dec_list = [-38.82364]

for ra, dec in zip(ra_list, dec_list):
    # ra, dec = 121.174886, -2.262535

    ticid, sector, cam, ccd = lcu.get_tess_lc(tess_dir, ra=ra, dec=dec)
    sector_dir = '/home/echickle/data/s%04d/s%04d-lc/'%(sector,sector)

    fname_atlas = lcu.get_atlas_lc(ticid, wd_tab, atlas_dir)
    fnames_ztf = lcu.get_ztf_lc(ra, dec)

    suffix = 'TIC_{}_ra_{}_dec_{}'.format(ticid, ra, dec)
    lcu.make_panel_plot(fname_atlas,fnames_ztf,sector_dir,ticid,cam,ccd,ra,dec,
                        gaia_tab,wd_tab,wd_main,rp_ext,out_dir,suffix,bls=bls,clip=clip)
