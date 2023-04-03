 # Download or locate ATLAS and ZTF light curves

import lc_utils as lcu
import os
import pdb

gaia_tab = "/home/echickle/data/100pc_clean.fits"
wd_tab = "/home/echickle/data/WDs.txt"
atlas_dir = "/data/ATLAS/"
out_dir = "/home/echickle/data/grnd_lc/"
out = open(out_dir + "s0061-cam3-ccd1.txt", "w")
out.write("TICID,CATALOG,FNAME\n")
fnames = ['wid_8_pow_9.779218_snr_0.25136_per_28.75983_TIC0000000768844169_cam_3_ccd_1_dur_0.025641026_epo_0.31623933_ra_115.83486156362_dec_-44.94315849534_phase_curve.png', 'wid_13_pow_9.54313_snr_0.19062_per_44.81933_TIC0000000820392808_cam_3_ccd_1_dur_0.041666668_epo_0.5694445_ra_126.70413453977_dec_-42.19810661004_phase_curve.png', 'wid_24_pow_9.490217_snr_0.10829_per_88.52323_TIC0000000822137226_cam_3_ccd_1_dur_0.14285715_epo_0.71428573_ra_124.74229868498_dec_-38.60304088786_phase_curve.png', 'wid_11_pow_9.407884_snr_0.19701_per_18.28732_TIC0000000819586980_cam_3_ccd_1_dur_0.033333335_epo_0.9333334_ra_131.80951922265_dec_-45.07048372993_phase_curve.png', 'wid_13_pow_9.438477_snr_0.15928_per_8.49752_TIC0000000822326111_cam_3_ccd_1_dur_0.055555556_epo_0.7037037_ra_121.59810917158_dec_-37.98794548443_phase_curve.png', 'wid_16_pow_11.38462_snr_0.18147_per_24.94436_TIC0000000813798045_cam_3_ccd_1_dur_0.05_epo_0.55_ra_127.89105486383_dec_-49.26028449771_phase_curve.png', 'wid_11_pow_9.971885_snr_0.13328_per_12.4869_TIC0000000822000394_cam_3_ccd_1_dur_0.07692308_epo_0.8974359_ra_123.3399345234_dec_-40.34215887177_phase_curve.png', 'wid_16_pow_9.819827_snr_0.16595_per_82.07464_TIC0000000822323963_cam_3_ccd_1_dur_0.06666667_epo_0.08888889_ra_121.39779581873_dec_-38.22592642046_phase_curve.png', 'wid_12_pow_9.991788_snr_0.16311_per_7.52789_TIC0000000768314250_cam_3_ccd_1_dur_0.033333335_epo_0.64444447_ra_117.66473254497_dec_-48.19590172737_phase_curve.png', 'wid_23_pow_11.210349_snr_0.13365_per_8.14466_TIC0000000822145894_cam_3_ccd_1_dur_0.1_epo_0.0_ra_124.23939677432_dec_-38.48430927063_phase_curve.png']

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
