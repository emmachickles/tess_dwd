import numpy as np
import os
import pandas as pd
import matplotlib.pyplot as plt
from astropy.coordinates import SkyCoord
import astropy.units as u 
from astroquery.gaia import Gaia
Gaia.ROW_LIMIT = 5
Gaia.MAIN_GAIA_TABLE = "gaiadr3.gaia_source"

mydir = '/home/echickle/data/EOF_IQR1010/'
vet_dir = mydir + 'vet_res/'
output_dir = mydir + 'met_plot/'
os.makedirs(output_dir, exist_ok=True)
gpu_dir = '/home/echickle/out/'
data_dir = '/home/echickle/data/'

# -- all -----------------------------------------------------------------------
pow_all, snr_all, wid_all, nt_all, dphi_all = [], [], [], [], []
gid_all, per_all = [], []
fnames = os.listdir(vet_dir)
for i in range(len(fnames)):
    # gid, pow, snr, wid, per, q, phi0, nt, dphi
    f = np.loadtxt(vet_dir + fnames[i])
    gid_all.extend(np.int64(f[:,0]))
    per_all.extend(f[:,4])
    pow_all.extend(f[:,1])
    snr_all.extend(f[:,2])
    wid_all.extend(np.int64(f[:,3]))
    nt_all.extend(np.int64(f[:,7]))
    dphi_all.extend(f[:,8])

fig1, ax1 = plt.subplots()
ax1.plot(pow_all, snr_all, '.k', ms=1, alpha=0.5)
ax1.set_xlabel('Peak Significance = (peak - median)/MAD')
ax1.set_ylabel('SNR = depth/MAD')
fig2, ax2 = plt.subplots()
ax2.plot(pow_all, wid_all, '.k', ms=1, alpha=0.5)
ax2.set_xlabel('Peak Significance = (peak - median)/MAD')
ax2.set_ylabel('Peak width')
fig3, ax3 = plt.subplots()
ax3.plot(nt_all, dphi_all, '.k', ms=1, alpha=0.5)
ax3.set_xlabel('Number of points in-eclipse')
ax3.set_ylabel('Max(diff(phi))')

wd_cat  = pd.read_csv(data_dir+'WDs.txt', header=None, sep='\s+', dtype='str')
ra = wd_cat[1].to_numpy().astype('float')
dec = wd_cat[2].to_numpy().astype('float')
gmag = wd_cat[4].to_numpy().astype('float')
source_id = wd_cat[3].to_numpy().astype('str')
gid_all, per_all = np.array(gid_all).astype('str'), np.array(per_all)
_, comm1, comm2 = np.intersect1d(gid_all, source_id, return_indices=True)
gmag_bkg, per_bkg = gmag[comm2], per_all[comm1]

plt.rcParams['agg.path.chunksize'] = 100000

# -- TESS ----------------------------------------------------------------------

f = np.loadtxt(gpu_dir+'GPU_TESS.result', usecols=(1,2,3,11,12), delimiter=',')
pow, snr, dphi = f[:,0], f[:,1], f[:,4]
wid, nt = np.int64(f[:,2]), np.int64(f[:,3])
match = np.loadtxt(gpu_dir+'GPU_TESS.result', usecols=(13), delimiter=',', dtype='str')
match = np.array([f[-1] for f in match])
inds = np.nonzero(match == 'T')
print('TESS Recovery: '+str(np.count_nonzero(match == 'T'))+' / '+str(len(match)))
ax1.plot(pow[inds], snr[inds], 'vr', label='TESS')
ax2.plot(pow[inds], wid[inds], 'vr', label='TESS')
ax3.plot(nt[inds], dphi[inds], 'vr', label='TESS')

# -- JVR -----------------------------------------------------------------------

# f cols: ra dec sig snr wid per q phi0 nt dphi
f = np.loadtxt(gpu_dir+'GPU_JVR.result', skiprows=1, usecols=(1,2,3,4,5,6,8,9,12,13))
ra, dec = f[:,0], f[:,1]
pow, snr, dphi = f[:,2], f[:,3], f[:,9]
wid, nt = np.int64(f[:,4]), np.int64(f[:,8])
per = f[:,5]

# cat cols : ra dec per
cat = np.loadtxt(data_dir+'ZTF_Eclipses.txt', usecols=(1,2,3))
co = SkyCoord(ra=ra, dec=dec, unit=(u.degree, u.degree), frame='icrs')
co_cat = SkyCoord(ra=cat[:,0], dec=cat[:,1], unit=(u.degree, u.degree), frame='icrs')
idx, d2d, _ = co.match_to_catalog_sky(co_cat)
good_idx = idx[np.nonzero(d2d.to(u.arcsec).value < 2)]

match = []
no_match = []
per_true = cat[:,2][good_idx]
dt = 2/1440
for i in range(len(per)):
    if (per[i] > per_true[i]-dt and per[i] < per_true[i]+dt) or\
    (per[i] > per_true[i]*2-dt and per[i] < per_true[i]*2+dt) or\
    (per[i] > per_true[i]/2-dt and per[i] < per_true[i]/2+dt):
        match.append(i)
    else:
        no_match.append(i)
match, no_match = np.array(match), np.array(no_match)

print('JVR Recovery: '+str(len(match))+' / '+str(len(per)))

ax1.plot(pow[match], snr[match], '<b', label='JVR')
ax2.plot(pow[match], wid[match], '<b', label='JVR')
ax3.plot(nt[match], dphi[match], '<b', label='JVR')

# gmag = []
# for i in range(len(co)):
#     j = Gaia.cone_search_async(co[i], radius=u.Quantity(3, u.arcsec))
#     if len(j.get_results()['phot_g_mean_mag']) > 0:
#         gmag.append(j.get_results()['phot_g_mean_mag'][0])
#     else:
#         gmag.append(np.nan)
# gmag = np.array(gmag)
# np.savetxt(gpu_dir+'Gmag_JVR.txt', gmag)
gmag = np.loadtxt(gpu_dir+'Gmag_JVR.txt')

fig4, ax4 = plt.subplots()
ax4.plot(gmag_bkg, per_bkg, '.k', ms=1,alpha=0.5)
ax4.set_xlabel('Gaia Gmag')
ax4.set_ylabel('Period [days]')
ax4.plot(gmag[match], per[match], '<b', label='JVR-recovered')
ax4.plot(gmag[no_match], per_true[no_match], 'Xr', label='JVR-unrecovered\nTrue period')
ax4.plot(gmag[no_match], per[no_match], 'Xm', label='JVR-unrecovered\nBLS period')
ax4.legend()
fig4.savefig(output_dir+'gmag_per_JVR.png', dpi=300)
ax4.set_ylim([0, 0.2])
fig4.savefig(output_dir+'gmag_per_JVR_zoom.png', dpi=300)

# -- KB UCB --------------------------------------------------------------------


# f cols: ra dec sig snr wid per q phi0 nt dphi
f = np.loadtxt(gpu_dir+'GPU_UCB.result', skiprows=1, usecols=(1,2,3,4,5,6,8,9,12,13))
ra, dec = f[:,0], f[:,1]
pow, snr, dphi = f[:,2], f[:,3], f[:,9]
wid, nt = np.int64(f[:,4]), np.int64(f[:,8])
per = f[:,5]
f = np.loadtxt(gpu_dir+'GPU_WDUCB.result', skiprows=1, usecols=(1,2,3,4,5,6,8,9,12,13))
ra = np.append(ra, f[:,0])
dec = np.append(dec, f[:,1])
pow = np.append(pow, f[:,2])
snr = np.append(snr, f[:,3])
dphi = np.append(dphi, f[:,9])
wid = np.append(wid, np.int64(f[:,4]))
nt = np.append(nt, np.int64(f[:,8]))
per = np.append(per, f[:,5])

cat = np.loadtxt(data_dir+"Kevin's UCBs - UCBs.csv", delimiter=',', skiprows=1, usecols=(2,3,4))
co = SkyCoord(ra=ra, dec=dec, unit=(u.degree, u.degree), frame='icrs')
co_cat = SkyCoord(ra=cat[:,0], dec=cat[:,1], unit=(u.degree, u.degree), frame='icrs')
idx, d2d, _ = co.match_to_catalog_sky(co_cat)
good_idx = idx[np.nonzero(d2d.to(u.arcsec).value < 2)]

match = []
no_match = []
per_true = cat[:,2][good_idx]
dt = 2/1440
for i in range(len(per)):
    if (per[i] > per_true[i]-dt and per[i] < per_true[i]+dt) or\
    (per[i] > per_true[i]*2-dt and per[i] < per_true[i]*2+dt) or\
    (per[i] > per_true[i]/2-dt and per[i] < per_true[i]/2+dt):
        match.append(i)
    else:
        no_match.append(i)
match, no_match = np.array(match), np.array(no_match)

print('KB UCB Recovery: '+str(len(match))+' / '+str(len(per)))

ax1.plot(pow[match], snr[match], '>g', label='KBUCB')
ax2.plot(pow[match], wid[match], '>g', label='KBUCB')
ax3.plot(nt[match], dphi[match], '>g', label='KBUCB')
# ax1.plot(pow[no_match], snr[no_match], 'Xg', label='KBUCB-unrecovered')
# ax2.plot(pow[no_match], wid[no_match], 'Xg', label='KBUCB-unrecovered')
# ax3.plot(nt[no_match], dphi[no_match], 'Xg', label='KBUCB-unrecovered')

# gmag = []
# for i in range(len(co)):
#     j = Gaia.cone_search_async(co[i], radius=u.Quantity(3, u.arcsec))
#     if len(j.get_results()['phot_g_mean_mag']) > 0:
#         gmag.append(j.get_results()['phot_g_mean_mag'][0])
#     else:
#         gmag.append(np.nan)
# gmag = np.array(gmag)
# np.savetxt(gpu_dir+'Gmag_KBUCB.txt', gmag)
gmag = np.loadtxt(gpu_dir+'Gmag_KBUCB.txt')

fig4, ax4 = plt.subplots()
ax4.plot(gmag_bkg, per_bkg, '.k', ms=1,alpha=0.5)
ax4.set_xlabel('Gaia Gmag')
ax4.set_ylabel('Period [days]')
ax4.plot(gmag[match], per[match], '>g', label='KBUCB-recovered')
ax4.plot(gmag[no_match], per_true[no_match], 'Xr', label='KBUCB-unrecovered\nTrue period')
ax4.plot(gmag[no_match], per[no_match], 'Xm', label='KBUCB-unrecovered\nBLS period')
ax4.legend()
fig4.savefig(output_dir+'gmag_per_KBUCB.png', dpi=300)
ax4.set_ylim([0, 0.1])
fig4.savefig(output_dir+'gmag_per_KBUCB_zoom.png', dpi=300)
       
# -- DWD -----------------------------------------------------------------------

# f cols: ra dec sig snr wid per q phi0 nt dphi
f = np.loadtxt(gpu_dir+'GPU_DWD.result', skiprows=1, usecols=(1,2,3,4,5,6,8,9,12,13))
ra, dec = f[:,0], f[:,1]
co = SkyCoord(ra=ra, dec=dec, unit=(u.degree, u.degree), frame='icrs')
pow, snr, dphi = f[:,2], f[:,3], f[:,9]
wid, nt = np.int64(f[:,4]), np.int64(f[:,8])
per = f[:,5]
match = np.array([0])
no_match = np.array([1,2,3,4])
print('Long DWD Recovery: '+str(len(match))+' / '+str(len(per)))
ax1.plot(pow[match], snr[match], 'oc', label='Long DWD')
ax2.plot(pow[match], wid[match], 'oc', label='Long DWD')
ax3.plot(nt[match], dphi[match], 'oc', label='Long DWD')

# gmag = []
# for i in range(len(co)):
#     try:
#         j = Gaia.cone_search_async(co[i], radius=u.Quantity(3, u.arcsec))
#     except:
#         gmag.append(np.nan)
#     if len(j.get_results()['phot_g_mean_mag']) > 0:
#         gmag.append(j.get_results()['phot_g_mean_mag'][0])
#     else:
#         gmag.append(np.nan)
# gmag = np.array(gmag)
# np.savetxt(gpu_dir+'Gmag_DWD.txt', gmag)
gmag = np.loadtxt(gpu_dir+'Gmag_DWD.txt')
per_true = np.array([0.11601549, 0.0800126, 0.2350606, 0.09986, 0.246137])

fig4, ax4 = plt.subplots()
ax4.plot(gmag_bkg, per_bkg, '.k', ms=1,alpha=0.5)
ax4.set_xlabel('Gaia Gmag')
ax4.set_ylabel('Period [days]')
ax4.plot(gmag[match], per_true[match], 'oc', label='Long DWD-recovered')
ax4.plot(gmag[no_match], per_true[no_match], 'Xr', label='Long DWD-unrecovered\nTrue period')
ax4.plot(gmag[no_match], per[no_match], 'Xm', label='Long DWD-unrecovered\nBLS period')
ax4.legend()
fig4.savefig(output_dir+'gmag_per_DWD.png', dpi=300)
ax4.set_ylim([0, 0.4])
fig4.savefig(output_dir+'gmag_per_DWD_zoom.png', dpi=300)

# -- save plots ----------------------------------------------------------------

ax1.legend()
fig1.tight_layout()
fig1.savefig(output_dir+'pow_snr.png', dpi=300)
print(output_dir+'pow_snr.png')
ax1.set_xlim([0,425])
ax1.set_ylim([0,15])
fig1.savefig(output_dir+'pow_snr_zoom.png', dpi=300)

ax2.legend()
fig2.tight_layout()
fig2.savefig(output_dir+'pow_wid.png', dpi=300)
print(output_dir+'pow_wid.png')
ax2.set_xlim([0,425])
ax2.set_ylim([0,175])
fig2.savefig(output_dir+'pow_wid_zoom.png', dpi=300)

ax3.legend()
fig3.tight_layout()
fig3.savefig(output_dir+'nt_dphi.png', dpi=300)
print(output_dir+'nt_dphi.png')   
ax3.set_xlim([0,600])
ax3.set_ylim([0,0.05])
fig3.savefig(output_dir+'nt_dphi_zoom.png', dpi=300)

