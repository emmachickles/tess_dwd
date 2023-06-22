import os
import numpy as np
import pandas as pd
from astropy.io import fits
import matplotlib.pyplot as plt
import lc_utils as lcu

sector = 60
data_dir = "/home/echickle/data/s%04d/"%sector+"s%04d-lc/"%sector
output_dir='/home/echickle/outlier_plots_TESS/'
os.makedirs(output_dir, exist_ok=True)
wd_tab = '/home/echickle/data/WDs.txt'

wd_cat = pd.read_csv(wd_tab, header=None, sep='\s+', dtype='str')
ra = wd_cat[1].to_numpy().astype('float')
dec = wd_cat[2].to_numpy().astype('float')
gmag = wd_cat[4].to_numpy().astype('float')
source_id = wd_cat[3].to_numpy().astype('str')
ticid = wd_cat[0].to_numpy().astype('int')

std_list = []
iqr_list = []
gmag_list = []
for cam in [1,2,3,4]:
    for ccd in [1,2,3,4]:
        suffix = "-{}-{}.npy".format(cam, ccd)
        y = np.load(data_dir+'lc'+suffix)
        t = np.load(data_dir+'ts'+suffix)
        ticid_ccd = np.load(data_dir+'id'+suffix).astype('int')
        _, _, comm2 = np.intersect1d(ticid_ccd, ticid, return_indices=True)
        gmag_list.extend(gmag[comm2])
        for i in range(len(y)):
            # >> remove nans
            inds = np.nonzero(~np.isnan(y[i]))
            t1, y1 = t[inds], y[i][inds]
            q3, q1 = np.percentile(y1, [75 ,25])
            iqr=(q3-q1)/2
            std_list.append(np.std(y1))
            iqr_list.append(iqr)

gmag_list, std_list, iqr_list = np.array(gmag_list), np.array(std_list), np.array(iqr_list)
np.savetxt(output_dir+'outlier_all.txt', np.array([gmag_list,iqr_list,std_list]))

# import astropy.units as u
# from astropy.coordinates import SkyCoord
# from astroquery.gaia import Gaia
# Gaia.ROW_LIMIT = 5
# Gaia.MAIN_GAIA_TABLE = "gaiadr3.gaia_source"

# data_dir = "/home/echickle/data/s%04d/"%sector+"s%04d-lc-ZTF/"%sector
# std_ztf = []
# iqr_ztf = []
# gmag_ztf = []
# for cam in [1,2,3,4]:
#     for ccd in [1,2,3,4]:
#         suffix = "-{}-{}.npy".format(cam, ccd)
#         print(suffix)
#         y = np.load(data_dir+'lc'+suffix)
#         t = np.load(data_dir+'ts'+suffix)
#         co_ccd = np.load(data_dir+'co'+suffix).astype('float')
        
#         for i in range(len(y)):
#             coord = SkyCoord(ra=co_ccd[i][0], dec=co_ccd[i][1],
#                              unit=(u.degree, u.degree), frame='icrs')
#             j = Gaia.cone_search_async(coord, radius=u.Quantity(3, u.arcsec))
#             if len(j.get_results()['phot_g_mean_mag']) > 0:
#                 gmag = j.get_results()['phot_g_mean_mag'][0]
#                 gmag_ztf.append(gmag)
            
#                 # >> remove nans
#                 inds = np.nonzero(~np.isnan(y[i]))
#                 t1, y1 = t[inds], y[i][inds]
#                 q3, q1 = np.percentile(y1, [75 ,25])
#                 iqr=(q3-q1)/2
#                 std_ztf.append(np.std(y1))
#                 iqr_ztf.append(iqr)

# gmag_ztf, std_ztf, iqr_ztf = np.array(gmag_ztf), np.array(std_ztf), np.array(iqr_ztf)
# np.savetxt(output_dir+'outlier_ztf.txt', np.array([gmag_ztf,iqr_ztf,std_ztf]))
data = np.loadtxt(output_dir + 'outlier_ztf.txt')
gmag_ztf, iqr_ztf, std_ztf = data[0], data[1], data[2]

plt.figure()
_ = plt.hist2d(gmag_list, std_list, bins=500, cmin=0.5)
plt.xlabel('Gmag')
plt.ylabel('Standard deviation')
plt.savefig(output_dir+'std_plot.png', dpi=300)
inds = np.nonzero(std_list < 50)
plt.figure()
_ = plt.hist2d(gmag_list[inds], std_list[inds], bins=500, cmin=0.5)
plt.plot(gmag_ztf, std_ztf, 'b<', label='JVR', ms=1)
plt.xlabel('Gmag')
plt.ylabel('Standard deviation')
plt.savefig(output_dir+'std_plot_zoom.png', dpi=300)

plt.figure()
_ = plt.hist2d(gmag_list, iqr_list, bins=500, cmin=0.5)
plt.xlabel('Gmag')
plt.ylabel('IQR')
plt.savefig(output_dir+'iqr_plot.png', dpi=300)
inds = np.nonzero(iqr_list < 25)
plt.figure()
_ = plt.hist2d(gmag_list[inds], iqr_list[inds], bins=500, cmin=0.5)
plt.plot(gmag_ztf, iqr_ztf, 'b<', label='JVR', ms=1)
plt.xlabel('Gmag')
plt.ylabel('IQR')
plt.savefig(output_dir+'iqr_plot_zoom.png', dpi=300)


