# Download or locate ATLAS and ZTF light curves

import lc_utils as lcu
import os

gaia_tab = "/home/echickle/data/100pc_clean.fits"
wd_tab = "/home/echickle/data/WDs.txt"
tess_dir = "/home/echickle/data/s0061/s0061-lc/"
atlas_dir = "/data/ATLAS/"
out_dir = "/home/echickle/data/grnd_lc/"
out = open(out_dir + "s0061-1-1.txt", "w")
out.write("TICID,CATALOG,FNAME\n")

fnames = ["wid_58_pow_72.73888_snr_0.53801_per_93.73281_TIC0000000800042858_cam_1_ccd_3_dur_0.1_epo_0.6_ra_134.442437099_dec_3.71537555667_phase_curve.png",
          "wid_37_pow_104.42451_snr_2.04829_per_138.58794_TIC0000000455206965_cam_1_ccd_4_dur_0.1_epo_0.16666667_ra_125.22306633609_dec_0.14541671625_phase_curve.png"]

ticid_list, ra_list, dec_list = [], [], []

for i in range(len(fnames)):
    fname = fnames[i].split("_")
    ticid, ra, dec = int(fname[8][3:]), float(fname[18]), float(fname[20])
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
