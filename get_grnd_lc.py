 # Download or locate ATLAS and ZTF light curves

import lc_utils as lcu
import os
import pdb

gaia_tab = "/home/echickle/data/100pc_clean.fits"
wd_tab = "/home/echickle/data/WDs.txt"
atlas_dir = "/data/ATLAS/"
out_dir = "/home/echickle/data/grnd_lc/"
os.makedirs(out_dir, exist_ok=True)

out = open(out_dir + "s0061-cam3-hsnr.txt", "w")
out.write("TICID,CATALOG,FNAME\n")

fnames = ['pow_29.4_snr_14.43_wid_68_per_64.082704_q_0.16666667_phi0_0.33333334_TIC0000000768424686_s0061_cam_3_ccd_1_ra_119.94515008557_dec_-46.54437568983__bls.png', 'pow_66.78_snr_17.17_wid_33_per_107.03394508_q_0.1_phi0_0.9666667_TIC0000000255756592_s0061_cam_3_ccd_3_ra_99.58194264646_dec_-48.98774139242__bls.png', 'pow_28.0_snr_10.94_wid_78_per_47.88138672_q_0.16666667_phi0_0.7777778_TIC0000000773409474_s0061_cam_3_ccd_2_ra_109.23865601627_dec_-37.47378523382__bls.png', 'pow_36.97_snr_15.64_wid_91_per_57.83196085_q_0.16666667_phi0_0.3888889_TIC0000000818613742_s0061_cam_3_ccd_1_ra_124.47274825261_dec_-45.1555678349__bls.png', 'pow_80.87_snr_18.46_wid_28_per_94.16322436_q_0.07692308_phi0_0.5641026_TIC0000000767706310_s0061_cam_3_ccd_3_ra_106.96774840196_dec_-47.43924221674__bls.png', 'pow_55.87_snr_26.68_wid_78_per_88.80873868_q_0.16666667_phi0_0.5555556_TIC0000000766421169_s0061_cam_3_ccd_4_ra_117.47092949719_dec_-55.61719603492__bls.png', 'pow_90.69_snr_28.46_wid_37_per_105.97455285_q_0.11111111_phi0_0.962963_TIC0000000096852226_s0061_cam_3_ccd_2_ra_105.07329557407_dec_-34.34172197073__bls.png', 'pow_57.57_snr_33.39_wid_56_per_185.50244233_q_0.16666667_phi0_0.5_TIC0000000811196429_s0061_cam_3_ccd_4_ra_126.74808455811_dec_-59.60437748422__bls.png']

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
