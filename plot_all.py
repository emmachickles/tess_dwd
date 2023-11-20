# ------------------------------------------------------------------------------

atlas_dir = "/data/ATLAS/WDs/"
out_dir = "/home/echickle/out/vet/"

# atlas_dir = "/data/ATLAS/sdB/"
# out_dir = "/home/echickle/out/vet_sdB/"

tess_dir = "/home/echickle/data/"
gaia_tab = "/home/echickle/data/100pc_clean.fits"
wd_tab = "/home/echickle/data/WDs.txt"
wd_main = "/data/GaiaEDR3_WD_main.fits"
rp_ext = "/data/GaiaEDR3_WD_RPM_ext.fits"

bls=True
clip=True

# ------------------------------------------------------------------------------

import lc_utils as lcu
import numpy as np
import sys
import os
import pdb
os.makedirs(out_dir, exist_ok=True)

# ------------------------------------------------------------------------------

# load data
ra_list = [36.4925618171062]
dec_list = [-69.3403824243844]
ra_list = [311.143166894]
dec_list = [-78.7005123396]
per_list = [None]*len(ra_list)

for ra, dec, per in zip(ra_list, dec_list, per_list):
    # ra, dec = 121.174886, -2.262535

    ticid, sector, cam, ccd = lcu.get_tess_lc(tess_dir, ra=ra, dec=dec)
    print('TIC '+str(ticid))
    print('sector '+str(sector))
    # sector,cam,ccd=56,4,4

    if ticid is None:
        sector_dir = None
    else:
        sector_dir = '/home/echickle/data/s%04d/s%04d-lc/'%(sector,sector)

    fname_atlas = lcu.get_atlas_lc(atlas_dir, ticid=ticid, wd_tab=wd_tab, ra=ra, dec=dec)
    fnames_ztf = lcu.get_ztf_lc(ra, dec)

    suffix = 'ra_{}_dec_{}'.format(ra, dec)
    lcu.make_panel_plot(fname_atlas,fnames_ztf,sector_dir,ticid,cam,ccd,ra,dec,
                        gaia_tab,wd_tab,wd_main,rp_ext,out_dir,suffix,bls=bls,clip=clip,
                        per=per, bins=70)
