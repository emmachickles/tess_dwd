import lc_utils as lcu
import sys


fname = "pow_72.52995_per_74.44577_TIC0000000803616096_cam_1_ccd_3_dur_0.01_epo_0.72999996_ra_131.26583606681_dec_-1.42705212618_phase_curve"
# fname = "pow_66.43712_per_25.71457_TIC0000000804368061_cam_1_ccd_3_dur_0.01_epo_0.39666665_ra_138.5436398888_dec_-1.21204428418_phase_curve"
# fname = sys.argv[1]
sector = 61

gaia_tab = "/home/echickle/data/100pc_clean.fits"
wd_tab = "/home/echickle/work/WDs.txt"
tess_dir = "/home/echickle/data/s0061/s0061-lc/"
atlas_dir = "/data/ATLAS/"
out_dir = "/home/echickle/data/s0061/s0061-plot/"

lcu.make_panel_plot(fname,sector,gaia_tab,wd_tab,tess_dir,atlas_dir,out_dir,bins=100)
