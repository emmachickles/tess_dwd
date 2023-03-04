import lc_utils as lcu
import sys
# fname = "pow_8.54768_per_14.46421_TIC0000000751298830_cam_1_ccd_1_dur_0.01_epo_0.9266666_ra_115.53650119612_dec_-7.56906251186_phase_curve.png"

# fname = sys.argv[1]
sector = 61

gaia_tab = "/home/echickle/data/100pc_clean.fits"
wd_tab = "/home/echickle/work/WDs.txt"
tess_dir = "/home/echickle/data/s0061/s0061-lc/"
atlas_dir = "/data/ATLAS/"
out_dir = "/home/echickle/data/s0061/s0061-plot/"

lcu.make_panel_plot(fname,sector,gaia_tab,wd_tab,tess_dir,atlas_dir,out_dir,bins=100)
