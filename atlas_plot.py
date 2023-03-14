import lc_utils as lcu
import sys
import os

sector = 61

gaia_tab = "/home/echickle/data/100pc_clean.fits"
wd_tab = "/home/echickle/work/WDs.txt"
tess_dir = "/home/echickle/data/s0061/s0061-lc/"
atlas_dir = "/data/ATLAS/"
out_dir = "/home/echickle/data/s0061/s0061-ls-1-2-230310-plot/" 
ls = True
os.makedirs(out_dir, exist_ok=True)

    # fname = "pow_66.43712_per_25.71457_TIC0000000804368061_cam_1_ccd_3_dur_0.01_epo_0.39666665_ra_138.5436398888_dec_-1.21204428418_phase_curve"
# fname = sys.argv[1]

fnames = ['pow_9.376388607962141_per_23.02528_TIC0000000836798156_cam_1_ccd_2_ra_127.76330633279_dec_-10.45665901568_phase_curve.png', 'pow_7.845325946450415_per_57.65564_TIC0000000836802291_cam_1_ccd_2_ra_127.55715940402_dec_-10.06556719683_phase_curve.png', 'pow_10.567752887232091_per_15.28579_TIC0000000837154537_cam_1_ccd_2_ra_129.56823793292_dec_-4.90766467546_phase_curve.png', 'pow_8.223340217472673_per_8.28835_TIC0000000835565870_cam_1_ccd_2_ra_127.74367875174_dec_-13.22022076108_phase_curve.png', 'pow_13.083769305640441_per_76.58879_TIC0000000836359453_cam_1_ccd_2_ra_133.08448001788_dec_-12.891784591_phase_curve.png', 'pow_11.54568101922127_per_7.95142_TIC0000000836683606_cam_1_ccd_2_ra_132.40797695836_dec_-10.73594152857_phase_curve.png', 'pow_7.795748617485594_per_50.38819_TIC0000000837042779_cam_1_ccd_2_ra_135.01794569597_dec_-6.70667118219_phase_curve.png', 'pow_14.150378464869366_per_73.04429_TIC0000000836628503_cam_1_ccd_2_ra_128.75152111111_dec_-11.52233740834_phase_curve.png', 'pow_6.689500817955809_per_38.7032_TIC0000000837027682_cam_1_ccd_2_ra_133.71656213501_dec_-7.56893737994_phase_curve.png', 'pow_10.609786782948788_per_15.82484_TIC0000000836509656_cam_1_ccd_2_ra_137.48211642696_dec_-9.70130113928_phase_curve.png', 'pow_10.886381712896796_per_12.05898_TIC0000000836932781_cam_1_ccd_2_ra_130.02368865862_dec_-7.52326284818_phase_curve.png', 'pow_10.859210353411399_per_9.33514_TIC0000000836719466_cam_1_ccd_2_ra_131.13975205423_dec_-10.00437637153_phase_curve.png']

for fname in fnames:

    lcu.make_panel_plot(fname,sector,gaia_tab,wd_tab,tess_dir,atlas_dir,out_dir,bins=100,ls=ls)

