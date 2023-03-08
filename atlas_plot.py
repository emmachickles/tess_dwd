import lc_utils as lcu
import sys
import os

sector = 61

gaia_tab = "/home/echickle/data/100pc_clean.fits"
wd_tab = "/home/echickle/work/WDs.txt"
tess_dir = "/home/echickle/data/s0061/s0061-lc/"
atlas_dir = "/data/ATLAS/"
# out_dir = "/home/echickle/data/s0061/s0061-plot/"
out_dir = "/home/echickle/data/s0061/s0061-bls-1-3-plot/" 
ls = False
os.makedirs(out_dir, exist_ok=True)

    # fname = "pow_66.43712_per_25.71457_TIC0000000804368061_cam_1_ccd_3_dur_0.01_epo_0.39666665_ra_138.5436398888_dec_-1.21204428418_phase_curve"
# fname = sys.argv[1]

fnames = ['pow_12.5425205_per_104.00709_TIC0000000842105238_cam_1_ccd_3_dur_0.022222223_epo_0.22222222_ra_141.24602051245_dec_0.29165606282_phase_curve.png', 'pow_24.177286_per_26.66217_TIC0000000800044379_cam_1_ccd_3_dur_0.01_epo_0.06666666_ra_135.459543098_dec_4.03645684226_phase_curve.png', 'pow_12.981104_per_38.33115_TIC0000000800145793_cam_1_ccd_3_dur_0.01_epo_0.55_ra_137.39155373749_dec_6.82227642745_phase_curve.png', 'pow_8.802645_per_21.29641_TIC0000000800073359_cam_1_ccd_3_dur_0.022222223_epo_0.26666668_ra_137.56217820053_dec_5.49667216936_phase_curve.png', 'pow_10.36024_per_17.86616_TIC0000000800073430_cam_1_ccd_3_dur_0.01_epo_0.37666667_ra_137.50715317601_dec_5.62623070292_phase_curve.png', 'pow_10.125452_per_14.41407_TIC0000000800288928_cam_1_ccd_3_dur_0.033333335_epo_0.10000001_ra_131.74917982566_dec_9.00812556033_phase_curve.png', 'pow_10.48829_per_37.64391_TIC0000000803812309_cam_1_ccd_3_dur_0.014925373_epo_0.8109453_ra_130.35263556567_dec_3.48029528955_phase_curve.png', 'pow_10.440945_per_13.45894_TIC0000000803781784_cam_1_ccd_3_dur_0.01_epo_0.6933333_ra_130.2336034432_dec_2.28816884181_phase_curve.png', 'pow_91.72858_per_93.73534_TIC0000000800042858_cam_1_ccd_3_dur_0.05_epo_0.51666665_ra_134.442437099_dec_3.71537555667_phase_curve.png', 'pow_28.573915_per_11.66773_TIC0000000804385658_cam_1_ccd_3_dur_0.014925373_epo_0.4925373_ra_136.4317652087_dec_-0.39988713417_phase_curve.png', 'pow_11.433591_per_13.68198_TIC0000000842214595_cam_1_ccd_3_dur_0.01_epo_0.5433333_ra_142.74030878624_dec_6.42160719161_phase_curve.png', 'pow_15.091732_per_26.66257_TIC0000000800180568_cam_1_ccd_3_dur_0.014925373_epo_0.11940298_ra_139.88515704792_dec_7.76020789802_phase_curve.png', 'pow_9.606545_per_8.47995_TIC0000000800174990_cam_1_ccd_3_dur_0.033333335_epo_0.7111111_ra_139.77132811702_dec_6.54630465641_phase_curve.png', 'pow_10.02875_per_10.83076_TIC0000000800066062_cam_1_ccd_3_dur_0.01_epo_0.22333333_ra_137.74469310714_dec_4.82624652883_phase_curve.png', 'pow_17.49422_per_39.33376_TIC0000000800108191_cam_1_ccd_3_dur_0.01_epo_0.20666666_ra_131.13347387363_dec_4.30275739454_phase_curve.png', 'pow_10.703911_per_83.56569_TIC0000000803619805_cam_1_ccd_3_dur_0.022222223_epo_0.43703705_ra_132.4267803417_dec_-0.82188395676_phase_curve.png', 'pow_10.46212_per_44.77879_TIC0000000842214509_cam_1_ccd_3_dur_0.014925373_epo_0.9004975_ra_142.63892913494_dec_6.31681550856_phase_curve.png', 'pow_9.606153_per_23.4388_TIC0000000803785483_cam_1_ccd_3_dur_0.01_epo_0.44_ra_131.18926172234_dec_2.81457187044_phase_curve.png', 'pow_10.570671_per_25.51211_TIC0000000842122942_cam_1_ccd_3_dur_0.022222223_epo_0.2_ra_140.58892372907_dec_2.27795467298_phase_curve.png']

for fname in fnames:

    lcu.make_panel_plot(fname,sector,gaia_tab,wd_tab,tess_dir,atlas_dir,out_dir,bins=100,ls=ls)
