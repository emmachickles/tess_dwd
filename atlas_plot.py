# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------

sector = 61
out_dir = "/home/echickle/data/s0061/s0061-bls-1-1-230322-plot/" 
ls = False

fnames = ['wid_29_pow_11.175058_snr_0.12944_per_25.64813_TIC0000000754957724_cam_1_ccd_1_dur_0.09090909_epo_0.51515156_ra_119.68515830855_dec_-3.04965004725_phase_curve.png', 'wid_14_pow_37.74368_snr_0.33903_per_8.88907_TIC0000000754984562_cam_1_ccd_1_dur_0.010638298_epo_0.8404255_ra_118.8871411275_dec_-1.78554163343_phase_curve.png', 'wid_15_pow_12.05825_snr_0.20198_per_9.5809_TIC0000000755007999_cam_1_ccd_1_dur_0.030303031_epo_0.3030303_ra_117.08183559093_dec_-2.67504248702_phase_curve.png', 'wid_19_pow_8.941142_snr_0.22476_per_24.71703_TIC0000000755224560_cam_1_ccd_1_dur_0.035714287_epo_0.23809525_ra_116.38898016943_dec_0.84949599802_phase_curve.png', 'wid_15_pow_8.76201_snr_0.16398_per_10.41067_TIC0000000755301147_cam_1_ccd_1_dur_0.05_epo_0.2_ra_116.32061061827_dec_0.95127080564_phase_curve.png']

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


for fname in fnames:

    lcu.make_panel_plot(fname,sector,gaia_tab,wd_tab,tess_dir,atlas_dir,out_dir,bins=100,ls=ls)

