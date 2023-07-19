# ------------------------------------------------------------------------------

atlas_dir = "/data/ATLAS/WDs/"
out_dir = "/home/echickle/out/vet/"

# atlas_dir = "/data/ATLAS/sdB/"
# out_dir = "/home/echickle/out/vet_sdB/"

tess_dir = "/home/echickle/data/"
gaia_tab = "/home/echickle/data/100pc_clean.fits"
wd_tab = "/home/echickle/data/WDs.txt"
wd_main = "/data/GaiaEDR3_WD_main.fits"
rp_ext = "/data/GaiaEDR3_WD_RPM_ext.fits"

bls=True
clip=True

# ------------------------------------------------------------------------------

import lc_utils as lcu
import numpy as np
import sys
import os
import pdb
os.makedirs(out_dir, exist_ok=True)

# ------------------------------------------------------------------------------

# load data
ra_list = [283.68673, 206.76688, 297.90905, 299.21942, 344.1389, 251.78027, 102.7166, 118.18445, 254.42919, 314.41071, 294.79548]
dec_list = [-35.51323, -31.14789, 34.12481, -31.65236, 47.70825, -24.91717, 53.12089, 64.39376, -19.40401, -46.54745, -40.57919]

ra_list = [344.1389, 314.41071, 297.90905, 254.42919, 251.78027, 206.76688, 118.18445]
dec_list = [47.70825, -46.54745, 34.12481, -19.40401, -24.91717, -31.14789, 64.39376]

# ra_list = [117.921569806686]
# dec_list = [-1.68915456479852]
# per_list =  [0.0800126]
# per = None
ra_list = [148.20492800852, 107.91352082077, 224.64672709097, 239.56511852792, 349.68735676108, 222.15113969775, 359.14113127532, 210.5430667278, 165.18471932946, 160.31838985217, 322.02920014796, 141.12244521273, 220.15073422962, 77.47048648895, 183.09137173963, 207.86173717535, 231.61248650513, 241.0831428595, 178.06537233829, 89.82399939144, 223.10909958156, 101.02979679798, 79.43643034424]
dec_list=[-38.13107858825, -32.04795874997, -51.09568874562, -64.96828301365, 47.68763795998, -49.7993870776, 41.98390420389, -67.76875503874, -58.07682245051, -41.38351485253, 33.08704661591, -44.15903759214, -65.73285566049, -69.51323232541, -49.27724379293, -30.92055972555, -41.92232392689, -32.80531640537, -42.98751233312, -63.75564377688, -70.73276488383, -54.85477838719, -69.6629598448]
per_list=[7.0212719, 9.71899294, 23.13373277, 19.59907167, 42.91262878, 18.70854436, 17.79024033, 7.14394663, 53.82189632, 14.82409192, 13.81254097, 46.89256953, 55.83774672, 12.56168812, 45.52471168, 6.76889178, 51.36198476, 13.71370881, 26.81552741, 38.15425364, 55.87896746, 42.0841303, 9.27063379]
per_list=np.array(per_list) / 1440.

ra_list=[121.1748862]
dec_list=[-2.262534927]
per_list=[None]

ra_list = [241.60995, 240.66688, 239.56512, 234.37716, 222.15114, 218.87704, 213.15951, 210.54307000000003, 174.54567, 160.31839, 145.76729, 122.7272]
dec_list = [-39.25712, -67.11267, -64.96828, -42.10545, -49.79939, -44.434059999999995, -48.60388, -67.76876, -51.66367, -41.38351, -31.95625, -47.48807]
per_list = [0.01478491, 0.01423518, 0.01361047, 0.01852404, 0.01299204, 0.02030136, 0.02021855, 0.00496107, 0.0096118, 0.01029451, 0.00831638, 0.01615914]

ra_list = [123.67249]
dec_list = [-64.44789]
per_list = [None]*len(ra_list)

for ra, dec, per in zip(ra_list, dec_list, per_list):
    # ra, dec = 121.174886, -2.262535

    ticid, sector, cam, ccd = lcu.get_tess_lc(tess_dir, ra=ra, dec=dec)
    print('TIC '+str(ticid))
    if ticid is None:
        sector_dir = None
    else:
        sector_dir = '/home/echickle/data/s%04d/s%04d-lc/'%(sector,sector)

    fname_atlas = lcu.get_atlas_lc(atlas_dir, ticid=ticid, wd_tab=wd_tab, ra=ra, dec=dec)
    fnames_ztf = lcu.get_ztf_lc(ra, dec)

    suffix = 'ra_{}_dec_{}'.format(ra, dec)
    lcu.make_panel_plot(fname_atlas,fnames_ztf,sector_dir,ticid,cam,ccd,ra,dec,
                        gaia_tab,wd_tab,wd_main,rp_ext,out_dir,suffix,bls=bls,clip=clip,
                        per=per, bins=70)
