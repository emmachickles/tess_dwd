import os
import numpy as np
import pandas as pd
from astropy.io import fits
import matplotlib.pyplot as plt
import lc_utils as lcu

data_dir='/data/ATLAS/'
skiprows=0
wd_tab = '/home/echickle/data/WDs.txt'
# wd_main='/data/GaiaEDR3_WD_main.fits'
# rp_ext='/data/GaiaEDR3_WD_RPM_ext.fits'

# source_id = np.empty(0)
# gmag = np.empty(0)

# maincat = fits.open(wd_main)
# source_id = np.append(source_id, maincat[1].data['source_id'])
# gmag = np.append(gmag, maincat[1].data['phot_g_mean_mag'])

# rpmext = fits.open(rp_ext)
# source_id = np.append(source_id, rpmext[1].data['source_id'])
# gmag = np.append(gmag, rpmext[1].data['phot_g_mean_mag'])

wd_cat = pd.read_csv(wd_tab, header=None, sep='\s+', dtype='str')
ra = wd_cat[1].to_numpy().astype('float')
dec = wd_cat[2].to_numpy().astype('float')
gmag = wd_cat[4].to_numpy().astype('float')
source_id = wd_cat[3].to_numpy().astype('str')

def process_lc(f):
    data=np.loadtxt(f,usecols=(0,3,4,16),skiprows=skiprows)
    Filter=np.loadtxt(f,usecols=(5),skiprows=skiprows,dtype=str)

    # >> filter by limiting magnitude                                           
    ZP=data[:,3] # >> zeropoint magnitude                                       
    Filter=Filter[ZP>17.5]
    data=data[ZP>17.5]
    t, y, dy = data[:,0], data[:,1], data[:,2]

    # >> remove nans                                                            
    inds = np.nonzero(~np.isnan(y))
    t, y, dy, Filter = t[inds], y[inds], dy[inds], Filter[inds]

    # >> match medians across filters (np.unique(Filter) = ['c', 'o'])          
    med = np.median(y[Filter == 'c'])
    med_o = np.median(y[Filter == 'o'])
    y[Filter == 'o'] += med - med_o

    # y, dy = lcu.normalize_lc(y, dy)
    
    q3, q1 = np.percentile(y, [75 ,25])
    iqr=(q3-q1)/2
    std = np.std(y)
    return iqr, std

# fnames = os.listdir(data_dir)
# id_list = []
# iqr_list = []
# std_list = []
# gmag_list = []
# for i in range(len(source_id)):
#     f = data_dir+source_id[i]    
#     if i % 1000 == 0:
#         print('{} / {}'.format(i, len(source_id)))

#     if os.path.exists(f):
#         iqr, std = process_lc(f)
#         id_list.append(np.int64(source_id[i]))
#         gmag_list.append(gmag[i])
#         iqr_list.append(iqr)
#         std_list.append(std)
# np.savetxt('/home/echickle/outlier.txt', np.array([id_list,gmag_list,iqr_list,std_list]))
# gmag_list, iqr_list = np.array(gmag_list), np.array(iqr_list)

data = np.loadtxt('/home/echickle/outlier.txt')
id_list, gmag_list, iqr_list, std_list = data[0], data[1], data[2], data[3]

# data = np.loadtxt('/home/echickle/iqr_unnorm.txt')
# gmag_list, iqr_list = data[:,0], data[:,1]


print(len(gmag_list))
from astropy.coordinates import SkyCoord
import astropy.units as u
data = np.loadtxt('/home/echickle/data/'+"Kevin's UCBs - UCBs.csv", delimiter=',', skiprows=1, usecols=(2,3,5))
ra_wd, dec_wd = data[:,0], data[:,1]
co = SkyCoord(ra=ra_wd, dec=dec_wd, unit=(u.degree, u.degree), frame='icrs')
cat = SkyCoord(ra=ra, dec=dec, unit=(u.degree, u.degree), frame='icrs')
idx, d2d, _ = co.match_to_catalog_sky(cat)
good_idx = idx[np.nonzero(d2d.to(u.arcsec).value < 2)]

gmag_wd, iqr_wd, std_wd = [], [], []
for i in range(len(good_idx)):
    f = data_dir+source_id[good_idx[i]]
    iqr, std = process_lc(f)
    gmag_wd.append(gmag[good_idx[i]])
    iqr_wd.append(iqr)
    std_wd.append(std)

    
cat_ztf = np.loadtxt("/home/echickle/data/ZTF_Eclipses.txt", dtype='str')
ra_ztf, dec_ztf = cat_ztf[:,1].astype('float'), cat_ztf[:,2].astype('float')
co = SkyCoord(ra=ra_ztf, dec=dec_ztf, unit=(u.degree, u.degree), frame='icrs')
cat = SkyCoord(ra=ra, dec=dec, unit=(u.degree, u.degree), frame='icrs')
idx, d2d, _ = co.match_to_catalog_sky(cat)
good_idx = idx[np.nonzero(d2d.to(u.arcsec).value < 2)]
gmag_ztf, iqr_ztf, std_ztf = [], [], []
for i in range(len(good_idx)):
    f = data_dir+source_id[good_idx[i]]
    
    if os.path.exists(f):
        iqr, std = process_lc(f)
        gmag_ztf.append(gmag[good_idx[i]])
        iqr_ztf.append(iqr)
        std_ztf.append(std)


# -- plotting ------------------------------------------------------------------
plt.figure()
_ = plt.hist2d(gmag_list, iqr_list, bins=1000, cmin=0.5)
plt.xlabel('Gmag')
plt.ylabel('IQR')
plt.savefig('/home/echickle/iqr_plot.png', dpi=300)
print('/home/echickle/iqr_plot.png')

inds = np.nonzero(iqr_list < 500)
plt.figure()
_ = plt.hist2d(gmag_list[inds], iqr_list[inds], bins=1000, cmin=0.5)
plt.xlabel('Gmag')
plt.ylabel('IQR')
plt.savefig('/home/echickle/iqr_plot_zoom.png', dpi=300)
print('/home/echickle/iqr_plot_zoom.png')

inds = np.nonzero(std_list < 500)
plt.figure()
_ = plt.hist2d(gmag_list[inds], std_list[inds], bins=1000, cmin=0.5)
plt.ylim([0,500])
plt.xlabel('Gmag')
plt.ylabel('STD')
# plt.plot(gmag_ztf, iqr_ztf, 'b<', label='JVR', ms=1)
# plt.plot(gmag_wd, iqr_wd, 'r>', label='UCB', ms=1)
plt.plot(gmag_ztf, std_ztf, 'b<', label='JVR', ms=1)
plt.plot(gmag_wd, std_wd, 'r>', label='UCB', ms=1)
plt.legend()
plt.savefig('/home/echickle/outlier_std.png', dpi=300)
