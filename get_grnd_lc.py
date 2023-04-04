 # Download or locate ATLAS and ZTF light curves

import lc_utils as lcu
import os
import pdb

gaia_tab = "/home/echickle/data/100pc_clean.fits"
wd_tab = "/home/echickle/data/WDs.txt"
atlas_dir = "/data/ATLAS/"
out_dir = "/home/echickle/data/grnd_lc/"
out = open(out_dir + "s0062-cam1.txt", "w")
out.write("TICID,CATALOG,FNAME\n")
fnames = ["pow_104.9595_snr_83.76016_wid_32_per_94.70348039_q_0.1_phi0_0.16666667_TIC0000000005393020_s0062_cam_1_ccd_1_ra_146.46253456217_dec_-19.73360846534_phase_curve.png",
          "pow_44.56791_snr_18.45331_wid_55_per_94.81745249_q_0.16666667_phi0_0.05555556_TIC0000000875661692_s0062_cam_1_ccd_1_ra_146.45323780121_dec_-19.73920357574_phase_curve.png",
          "pow_28.372883_snr_10.10433_wid_11_per_85.12449224_q_0.03846154_phi0_0.7820513_TIC0000000875850017_s0062_cam_1_ccd_2_ra_150.3993200987_dec_-17.65750349749_phase_curve.png",
          "pow_63.975212_snr_19.23714_wid_63_per_167.57943106_q_0.16666667_phi0_0.2777778_TIC0000000471013677_s0062_cam_1_ccd_3_ra_156.457125_dec_0.65170833_phase_curve.png",
          "pow_78.820274_snr_17.78026_wid_21_per_148.82253029_q_0.055555556_phi0_0.2777778_TIC0000000101433897_s0062_cam_1_ccd_2_ra_149.05450015443_dec_-19.51482476495_phase_curve.png",
          "pow_28.372883_snr_10.10433_wid_11_per_85.12449224_q_0.03846154_phi0_0.7820513_TIC0000000875850017_s0062_cam_1_ccd_2_ra_150.3993200987_dec_-17.65750349749_phase_curve.png",
          "pow_38.45254_snr_13.03353_wid_60_per_99.97842552_q_0.16666667_phi0_0.6111111_TIC0000000875790341_s0062_cam_1_ccd_2_ra_150.54888973454_dec_-19.4270185001_phase_curve.png",
          "pow_74.32736_snr_24.1244_wid_39_per_133.89990271_q_0.06666667_phi0_0.08888889_TIC0000000471014834_s0062_cam_1_ccd_3_ra_156.722872_dec_-10.225034_phase_curve.png"]

for i in range(len(fnames)):
    fname = fnames[i].split("_")
    ticid, ra, dec = int(fname[12][3:]), float(fname[19]), float(fname[21])
    fname = "_".join(fname)
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
