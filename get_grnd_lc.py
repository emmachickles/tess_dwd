 # Download or locate ATLAS and ZTF light curves

import lc_utils as lcu
import os
import pdb

gaia_tab = "/home/echickle/data/100pc_clean.fits"
wd_tab = "/home/echickle/data/WDs.txt"
atlas_dir = "/data/ATLAS/"
out_dir = "/home/echickle/data/mk_lc/"
os.makedirs(out_dir, exist_ok=True)

out = open(out_dir + "mk-list.txt", "w") # !!
out.write("TICID,CATALOG,FNAME\n")

# fnames = ['pow_108.88_snr_98.63_wid_41_per_170.64528885_q_0.082_phi0_0.49614_TIC0000000193092806_s0062_cam_2_ccd_3_ra_156.73531121977_dec_-27.38252733725__bls.png', 'pow_57.7_snr_14.42_wid_31_per_126.62865026_q_0.108_phi0_0.07966_TIC0000000808364853_s0062_cam_3_ccd_4_ra_129.72499594618_dec_-59.72792439214__bls.png', 'pow_29.1_snr_9.82_wid_11_per_85.12449224_q_0.05_phi0_0.77844_TIC0000000875850017_s0062_cam_1_ccd_2_ra_150.3993200987_dec_-17.65750349749__bls.png', 'pow_33.3_snr_11.06_wid_24_per_108.65613705_q_0.055_phi0_0.84401_TIC0000000826164830_s0062_cam_2_ccd_1_ra_133.57593122817_dec_-40.03515344011__bls.png', 'pow_79.64_snr_18.55_wid_29_per_94.1608029_q_0.087_phi0_0.55265_TIC0000000767706310_s0061_cam_3_ccd_3_ra_106.96774840196_dec_-47.43924221674__bls.png']

ticid_list = [767706310, 675125497, 808364853]
ra_list = [106.967748, 72.3806833, 129.7249958]
dec_list = [-47.439242, -78.4414972, -59.7279250]
for i in range(len(ticid_list)):

# for i in range(len(fnames)):
    # fname = fnames[i].split("_")
    # ticid, ra, dec = int(fname[12][3:]), float(fname[19]), float(fname[21])
    # fname = "_".join(fname)

    fname = ""
    ticid, ra, dec = ticid_list[i], ra_list[i], dec_list[i]
    
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
