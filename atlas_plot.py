# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------

sector = 61
out_dir = "/home/echickle/data/s0061/s0061-bls-1-1-230322-plot/" 
ls = False

fnames = ['wid_29_pow_11.175058_snr_0.12944_per_25.64813_TIC0000000754957724_cam_1_ccd_1_dur_0.09090909_epo_0.51515156_ra_119.68515830855_dec_-3.04965004725_phase_curve.png', 'wid_12_pow_8.686864_snr_0.23057_per_70.85747_TIC0000000755200458_cam_1_ccd_1_dur_0.02173913_epo_0.02173913_ra_117.51599517275_dec_0.12763613213_phase_curve.png', 'wid_9_pow_9.465144_snr_0.25537_per_10.7292_TIC0000000755024789_cam_1_ccd_1_dur_0.023809524_epo_0.42857143_ra_115.63168585913_dec_-2.22510952883_phase_curve.png', 'wid_11_pow_9.574354_snr_0.20893_per_12.35618_TIC0000000755054313_cam_1_ccd_1_dur_0.05263158_epo_0.22807018_ra_118.79262186423_dec_-1.5264232671_phase_curve.png', 'wid_10_pow_117.502975_snr_0.90613_per_174.2141_TIC0000000803489769_cam_1_ccd_1_dur_0.023809524_epo_0.21428572_ra_121.17488619979_dec_-2.26253492655_phase_curve.png', 'wid_14_pow_37.74368_snr_0.33903_per_8.88907_TIC0000000754984562_cam_1_ccd_1_dur_0.010638298_epo_0.8404255_ra_118.8871411275_dec_-1.78554163343_phase_curve.png', 'wid_15_pow_12.05825_snr_0.20198_per_9.5809_TIC0000000755007999_cam_1_ccd_1_dur_0.030303031_epo_0.3030303_ra_117.08183559093_dec_-2.67504248702_phase_curve.png', 'wid_11_pow_9.940201_snr_0.17192_per_10.86086_TIC0000000836033167_cam_1_ccd_1_dur_0.045454547_epo_0.45454547_ra_123.38667416993_dec_-11.4702153666_phase_curve.png', 'wid_10_pow_10.829989_snr_0.27439_per_9.07523_TIC0000000755214219_cam_1_ccd_1_dur_0.025641026_epo_0.7692308_ra_116.19981624721_dec_0.36631943402_phase_curve.png', 'wid_8_pow_8.88976_snr_0.24563_per_7.44588_TIC0000000803852871_cam_1_ccd_1_dur_0.02173913_epo_0.7173913_ra_120.04952985372_dec_-1.45654140585_phase_curve.png', 'wid_19_pow_8.941142_snr_0.22476_per_24.71703_TIC0000000755224560_cam_1_ccd_1_dur_0.035714287_epo_0.23809525_ra_116.38898016943_dec_0.84949599802_phase_curve.png', 'wid_10_pow_9.223532_snr_0.16791_per_7.71867_TIC0000000755035132_cam_1_ccd_1_dur_0.035714287_epo_0.1904762_ra_117.05053463845_dec_-1.52636694638_phase_curve.png', 'wid_12_pow_9.52494_snr_0.20287_per_7.62223_TIC0000000755216985_cam_1_ccd_1_dur_0.05882353_epo_0.21568628_ra_116.99996472621_dec_0.31741312848_phase_curve.png', 'wid_15_pow_8.76201_snr_0.16398_per_10.41067_TIC0000000755301147_cam_1_ccd_1_dur_0.05_epo_0.2_ra_116.32061061827_dec_0.95127080564_phase_curve.png', 'wid_8_pow_10.825268_snr_0.24253_per_8.19767_TIC0000000755056336_cam_1_ccd_1_dur_0.027777778_epo_0.37962964_ra_119.22734638087_dec_-1.42351689182_phase_curve.png', 'wid_11_pow_8.507903_snr_0.19457_per_47.61053_TIC0000000803454388_cam_1_ccd_1_dur_0.02173913_epo_0.6304348_ra_126.11537307001_dec_-2.01659526362_phase_curve.png']

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

