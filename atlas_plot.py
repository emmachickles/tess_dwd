import lc_utils as lcu
import sys
import os

sector = 57

gaia_tab = "/home/echickle/data/100pc_clean.fits"
wd_tab = "/home/echickle/work/WDs.txt"
tess_dir = "/home/echickle/data/s0057/s0057-lc/"
atlas_dir = "/data/ATLAS/"
# out_dir = "/home/echickle/data/s0061/s0061-plot/"
out_dir = "/home/echickle/data/s0057/s0057-bls-4-3-plot/" 
ls = False
os.makedirs(out_dir, exist_ok=True)

    # fname = "pow_66.43712_per_25.71457_TIC0000000804368061_cam_1_ccd_3_dur_0.01_epo_0.39666665_ra_138.5436398888_dec_-1.21204428418_phase_curve"
# fname = sys.argv[1]

fnames = ['pow_10.388211_per_14.39164_TIC0000001102463585_cam_4_ccd_3_dur_0.014925373_epo_0.840796_ra_233.36574765131_dec_65.63299107845_phase_curve.png', 'pow_8.6986685_per_9.87864_TIC0000001102562633_cam_4_ccd_3_dur_0.014925373_epo_0.48258704_ra_235.02460635137_dec_71.42844668605_phase_curve.png', 'pow_10.258971_per_8.24602_TIC0000001271224319_cam_4_ccd_3_dur_0.01_epo_0.7633333_ra_256.04661400869_dec_64.50173472853_phase_curve.png', 'pow_8.710852_per_17.70584_TIC0000001271299893_cam_4_ccd_3_dur_0.01_epo_0.15666667_ra_252.07241837187_dec_71.17570220825_phase_curve.png', 'pow_8.903435_per_22.64025_TIC0000001271324893_cam_4_ccd_3_dur_0.014925373_epo_0.840796_ra_251.06433982368_dec_73.21530511774_phase_curve.png', 'pow_10.135358_per_7.19228_TIC0000001102469580_cam_4_ccd_3_dur_0.022222223_epo_0.0074074077_ra_234.44654357109_dec_67.32423744556_phase_curve.png', 'pow_10.519142_per_13.85007_TIC0000001201351162_cam_4_ccd_3_dur_0.033333335_epo_0.7444445_ra_240.47840328564_dec_73.15150887754_phase_curve.png', 'pow_10.140496_per_61.36082_TIC0000000232605099_cam_4_ccd_3_dur_0.05_epo_0.016666668_ra_257.92671356738_dec_64.47981844148_phase_curve.png', 'pow_8.700124_per_16.15434_TIC0000001201215001_cam_4_ccd_3_dur_0.01_epo_0.11_ra_249.33910825166_dec_62.25983247138_phase_curve.png', 'pow_10.233966_per_78.94593_TIC0000001201270515_cam_4_ccd_3_dur_0.01_epo_0.48666665_ra_243.45764631037_dec_65.00434638116_phase_curve.png', 'pow_10.366046_per_33.33109_TIC0000001102441451_cam_4_ccd_3_dur_0.033333335_epo_0.92222226_ra_237.41058193869_dec_64.62746569031_phase_curve.png', 'pow_10.212882_per_7.01903_TIC0000001102495641_cam_4_ccd_3_dur_0.01_epo_0.8666667_ra_236.73675008043_dec_68.89687512805_phase_curve.png', 'pow_10.195262_per_45.85868_TIC0000001201347741_cam_4_ccd_3_dur_0.01_epo_0.37666667_ra_249.03747632316_dec_73.95950826764_phase_curve.png', 'pow_8.70013_per_9.06975_TIC0000001102427704_cam_4_ccd_3_dur_0.014925373_epo_0.6069652_ra_236.15895409038_dec_63.10121497308_phase_curve.png', 'pow_9.33968_per_11.24652_TIC0000001271298334_cam_4_ccd_3_dur_0.014925373_epo_0.85572135_ra_255.4988279216_dec_71.37412823686_phase_curve.png', 'pow_10.188189_per_9.85554_TIC0000001271209436_cam_4_ccd_3_dur_0.01_epo_0.42_ra_257.9615941421_dec_63.81627789551_phase_curve.png']

for fname in fnames:

    lcu.make_panel_plot(fname,sector,gaia_tab,wd_tab,tess_dir,atlas_dir,out_dir,bins=100,ls=ls)

