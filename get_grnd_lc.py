 # Download or locate ATLAS and ZTF light curves

import lc_utils as lcu
import os
import pdb

gaia_tab = "/home/echickle/data/100pc_clean.fits"
wd_tab = "/home/echickle/data/WDs.txt"
atlas_dir = "/data/ATLAS/"
out_dir = "/home/echickle/data/grnd_lc/"
os.makedirs(out_dir, exist_ok=True)

out = open(out_dir + "s0062-cam3-hsnr.txt", "w")
out.write("TICID,CATALOG,FNAME\n")

fnames = ['pow_22.21_snr_10.15_wid_109_per_136.15979807_q_0.16666667_phi0_0.33333334_TIC0000000812106363_s0062_cam_3_ccd_1_ra_133.01229323239_dec_-54.38831841839__bls.png', 'pow_57.86_snr_27.57_wid_72_per_64.08902868_q_0.16666667_phi0_0.8333334_TIC0000000768424686_s0062_cam_3_ccd_2_ra_119.94515008557_dec_-46.54437568983__bls.png', 'pow_19.87_snr_10.06_wid_58_per_28.79336164_q_0.16666667_phi0_0.9444445_TIC0000000817353844_s0062_cam_3_ccd_2_ra_121.36001730254_dec_-48.75439361901__bls.png', 'pow_42.72_snr_16.69_wid_59_per_93.52429045_q_0.16666667_phi0_0.6666667_TIC0000000858457131_s0062_cam_3_ccd_1_ra_145.75785809212_dec_-54.47926819596__bls.png', 'pow_28.06_snr_12.5_wid_84_per_142.6865767_q_0.16666667_phi0_0.72222227_TIC0000000808487411_s0062_cam_3_ccd_4_ra_130.04401743659_dec_-58.02920859833__bls.png', 'pow_30.63_snr_19.91_wid_279_per_98.74142425_q_0.16666667_phi0_0.6666667_TIC0000000808091787_s0062_cam_3_ccd_4_ra_132.80334515325_dec_-61.7585621099__bls.png', 'pow_28.51_snr_11.77_wid_75_per_77.5724891_q_0.16666667_phi0_0.22222222_TIC0000000806495373_s0062_cam_3_ccd_4_ra_125.41921350374_dec_-61.94281657399__bls.png', 'pow_59.7_snr_36.0_wid_56_per_185.51190601_q_0.16666667_phi0_0.3888889_TIC0000000811196429_s0062_cam_3_ccd_4_ra_126.74808455811_dec_-59.60437748422__bls.png', 'pow_31.78_snr_11.49_wid_66_per_114.55179436_q_0.16666667_phi0_0.16666667_TIC0000000816100516_s0062_cam_3_ccd_1_ra_139.61803463465_dec_-45.61971325733__bls.png', 'pow_19.62_snr_10.04_wid_119_per_72.58792115_q_0.16666667_phi0_0.8333334_TIC0000000810710918_s0062_cam_3_ccd_1_ra_138.72312506737_dec_-54.30828732314__bls.png', 'pow_58.12_snr_14.69_wid_28_per_126.63298064_q_0.071428575_phi0_0.0952381_TIC0000000808364853_s0062_cam_3_ccd_4_ra_129.72499594618_dec_-59.72792439214__bls.png', 'pow_23.49_snr_13.34_wid_74_per_83.44225727_q_0.16666667_phi0_0.5555556_TIC0000000859096828_s0062_cam_3_ccd_1_ra_140.90679073595_dec_-54.60157022607__bls.png', 'pow_20.23_snr_10.07_wid_132_per_94.40409889_q_0.16666667_phi0_0.9444445_TIC0000000869413651_s0062_cam_3_ccd_1_ra_144.92906868607_dec_-45.73147179484__bls.png']

# ticid_list = [767706310, 803489769]
# ra_list = [106.967748, 121.174886]
# dec_list = [-47.439242, -2.262535]
# for i in range(len(ticid_list)):

for i in range(len(fnames)):
    fname = fnames[i].split("_")
    ticid, ra, dec = int(fname[12][3:]), float(fname[19]), float(fname[21])
    fname = "_".join(fname)

    # fname = ""
    # ticid, ra, dec = ticid_list[i], ra_list[i], dec_list[i]
    
    out.write("{},{},{}\n".format(ticid, "TESS", fname))
    
    fpath = lcu.get_atlas_lc(ticid, wd_tab, atlas_dir)
    if type(fpath) == type(""):
        fname = fpath.split("/")[-1]
        os.system("cp "+fpath+" "+out_dir+fname)
        out.write("{},{},{}\n".format(ticid, "ATLAS", fname))

    fpaths = lcu.get_ztf_lc(ra, dec)
    for j in range(len(fpaths)):
        fpath = fpaths[j]
        fname = fpath.split("/")[-1]
        os.system("mv "+fpath+" "+out_dir+fname)
        out.write("{},{},{}\n".format(ticid, "ZTF", fname))    

out.close()
