# ------------------------------------------------------------------------------

# atlas_dir = "/data/ATLAS/WDs/"
# out_dir = "/home/echickle/out/vet/"

atlas_dir = "/data/ATLAS/sdB/"
out_dir = "/home/echickle/out/vet_sdB/"

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
ra_list = [283.68673, 206.76688, 297.90905, 299.21942, 344.1389, 251.78027, 102.7166, 118.18445, 254.42919, 314.41071, 294.79548]
dec_list = [-35.51323, -31.14789, 34.12481, -31.65236, 47.70825, -24.91717, 53.12089, 64.39376, -19.40401, -46.54745, -40.57919]

ra_list = [118.18445, 254.42919, 314.41071, 294.79548]
dec_list =[64.39376, -19.40401, -46.54745, -40.57919]


for ra, dec in zip(ra_list, dec_list):
    # ra, dec = 121.174886, -2.262535

    ticid, sector, cam, ccd = lcu.get_tess_lc(tess_dir, ra=ra, dec=dec)
    if ticid is None:
        sector_dir = None
    else:
        sector_dir = '/home/echickle/data/s%04d/s%04d-lc/'%(sector,sector)

    fname_atlas = lcu.get_atlas_lc(atlas_dir, ticid=ticid, wd_tab=wd_tab, ra=ra, dec=dec)
    fnames_ztf = lcu.get_ztf_lc(ra, dec)

    suffix = 'ra_{}_dec_{}'.format(ra, dec)
    lcu.make_panel_plot(fname_atlas,fnames_ztf,sector_dir,ticid,cam,ccd,ra,dec,
                        gaia_tab,wd_tab,wd_main,rp_ext,out_dir,suffix,bls=bls,clip=clip)
