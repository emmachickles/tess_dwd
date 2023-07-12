import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from astropy.coordinates import SkyCoord
import astropy.units as u 
from astroquery.gaia import Gaia
import lc_utils as lcu
from vet_utils import *
import pdb

# Set constants 
Gaia.ROW_LIMIT = 5
Gaia.MAIN_GAIA_TABLE = 'gaiadr3.gaia_source'
out_dir = '/scratch/echickle/tess/vet/'
data_dir = '/scratch/echickle/tess/BLS_results/'
lc_dir = '/scratch/data/tess/lcur/ffi/'
qflag_dir = '/scratch/echickle/QLPqflags/'
wd_tab = '/scratch/echickle/WDs.txt'
wd_main = '/scratch/echickle/GaiaEDR3_WD_main.fits'
rp_ext = '/scratch/echickle/GaiaEDR3_WD_RPM_ext.fits'
sector_list = [56,57,58,59,60,61,62,63,64,65]
pmax = 0.15*1440

# Make output directory
os.makedirs(out_dir, exist_ok=True)

# Load Gaia white dwarf results 
result_list = append_result_file(data_dir, sector_list)

# Initialize metric plots 
fig1, ax1 = plt.subplots(figsize=(5,4)) # Power vs. SNR
ax1.plot(result_list[:,3], result_list[:,4], '.k', ms=1, alpha=0.5)

fig2, ax2 = plt.subplots(figsize=(5,4)) # Power vs. width
ax2.plot(result_list[:,3], result_list[:,5], '.k', ms=1, alpha=0.5)

fig3, ax3 = plt.subplots(figsize=(5,4)) # nTransit vs dphi
ax3.plot(result_list[:,7], result_list[:,8], '.k', ms=1, alpha=0.5)

# -- catalogs -----------------------------------------------------------------

# >> signals discovered in S61
ticid_tess = np.array([803489769, 36085812, 800042858, 270679102, 455206965, 452954413, 767706310, 96852226, 677736827])
ticid_tess, result_tess = match_catalog(result_list, ticid_tess)
ax1.plot(result_tess[:,3], result_tess[:,4], 'vr', label='TESS')
ax2.plot(result_tess[:,3], result_tess[:,5], 'vr', label='TESS')
ax3.plot(result_tess[:,7], result_tess[:,8], 'vr', label='TESS')

# >> KB UCBs systems in the white dwarf catalog
ticid_ucb = np.array([719830487, 702774311, 713144139, 1939638541, 746047312, 1504732386, 1717410027, 1958881490, 2040677137, 1270451609])
period_ucb_true = np.array([0.010030088, 0.014240466, 0.014282094, 0.0144914, 0.016464692, 0.018356874, 0.0281957, 0.029969971, 0.038365738, 0.04542631929])
ticid_ucb, period_ucb_true, result_ucb = match_catalog(result_list, ticid_ucb, period_ucb_true)
match_ucb = match_period(result_ucb[:,6], period_ucb_true)
if np.count_nonzero(match_ucb) > 0:
    ax1.plot(result_ucb[:,3][match_ucb], result_ucb[:,4][match_ucb], '>g', label='UCB')
    ax2.plot(result_ucb[:,3][match_ucb], result_ucb[:,5][match_ucb], '>g', label='UCB')
    ax3.plot(result_ucb[:,7][match_ucb], result_ucb[:,8][match_ucb], '>g', label='UCB')

# >> long period DWDs in the white dwarf catalog
ticid_dwd = np.array([755049150, 903242599, 840220096, 219868627])
period_dwd_true = np.array([0.0800126, 0.09986, 0.11601549, 0.246137])
ticid_dwd, period_dwd_true, result_dwd = match_catalog(result_list, ticid_dwd, period_dwd_true)
match_dwd = match_period(result_dwd[:,6], period_dwd_true)

if np.count_nonzero(match_dwd) > 0:
    ax1.plot(result_dwd[:,3][match_dwd], result_dwd[:,4][match_dwd], 'oc', label='DWD')
    ax2.plot(result_dwd[:,3][match_dwd], result_dwd[:,5][match_dwd], 'oc', label='DWD')
    ax3.plot(result_dwd[:,7][match_dwd], result_dwd[:,8][match_dwd], 'oc', label='DWD')

print('KB UCB Recovery: '+str(np.count_nonzero(match_ucb))+' / '+str(len(match_ucb)))
print('Recovered UCBs: '+','.join(ticid_ucb[match_ucb].astype('str')))
print('DWD Recovery: '+str(np.count_nonzero(match_dwd))+' / '+str(len(match_dwd)))
print('Recovered DWDs: '+','.join(ticid_dwd[match_dwd].astype('str')))

# -- Gmag plots ---------------------------------------------------------------

plot_gmag(out_dir, wd_tab, result_list, result_ucb, match_ucb, period_ucb_true, suffix='UCB')
plot_gmag(out_dir, wd_tab, result_list, result_dwd, match_dwd, period_dwd_true, suffix='DWD')

# -- JVR -----------------------------------------------------------------------

data_dir = '/scratch/echickle/tess_jvr/BLS_results/'

# Load Jan's WDRD catalog
cat = np.loadtxt("/scratch/echickle/ZTF_Eclipses.txt", usecols=(1,2,3))
ra_ztf, dec_ztf, period_ztf_true = cat[:,0], cat[:,1], cat[:,2]

# Load Jan's WDRD results
result_ztf = append_result_file(data_dir, sector_list)

ra_ztf, dec_ztf, period_ztf_true, result_ztf = match_coord(ra_ztf, dec_ztf, period_ztf_true, result_ztf)
match_ztf = match_period(result_ztf[:,6], period_ztf_true)
if np.count_nonzero(match_ztf) > 0:
    ax1.plot(result_ztf[:,3][match_ztf], result_ztf[:,4][match_ztf], '<b', label='JVR')
    ax2.plot(result_ztf[:,3][match_ztf], result_ztf[:,5][match_ztf], '<b', label='JVR')
    ax3.plot(result_ztf[:,7][match_ztf], result_ztf[:,8][match_ztf], '<b', label='JVR')

print('DWD Recovery: '+str(np.count_nonzero(match_ztf))+' / '+str(len(match_ztf)))

gmag_ztf = get_gmag(result_ztf, data_dir)
plot_gmag(out_dir, wd_tab, result_list, result_ztf, match_ztf, period_ztf_true,
          gmag_catalog=gmag_ztf, suffix='JVR')

# -- save figures --------------------------------------------------------------

ax1.set_xlabel('Peak Significance = (peak - median)/MAD')
ax1.set_ylabel('SNR = depth/MAD')
ax1.legend()
fig1.tight_layout()
fig1.savefig(out_dir + 'pow_snr.png', dpi=300)
print(out_dir + 'pow_snr.png')
ax1.set_xlim([0, 7500])
ax1.set_ylim([0, 500])
fig1.savefig(out_dir + 'pow_snr_zoom.png', dpi=300)
print(out_dir + 'pow_snr_zoom.png')

ax2.set_xlabel('Peak Significance = (peak - median)/MAD')
ax2.set_ylabel('Peak width')
ax2.legend()
fig2.tight_layout()
fig2.savefig(out_dir + 'pow_wid.png', dpi=300)
print(out_dir + 'pow_wid.png')
ax2.set_xlim([0, 7500])
ax2.set_ylim([0, 100])
fig2.savefig(out_dir + 'pow_wid_zoom.png', dpi=300)
print(out_dir + 'pow_wid_zoom.png')

ax3.set_xlabel('Number of points in transit')
ax3.set_ylabel('Dphi')
ax3.legend()
fig3.tight_layout()
fig3.savefig(out_dir + 'nt_dphi.png', dpi=300)
print(out_dir + 'nt_dphi.png')

ax3.set_ylim([-0.005, 0.04])
fig3.savefig(out_dir + 'nt_dphi_zoom.png', dpi=300)

# ------------------------------------------------------------------------------

# ra, dec, pow, snr, dphi, wid, nt, per = [np.empty(0)]*8
# fnames_jvr = []

# # ticid, ra, dec, sig, snr, wid, period, period_min, q, phi0, epo, rp, nt, dphi
# for sector in [56]: # !!
#     res_dir = data_dir + 's%04d-ZTF/'%sector + 's%04d-gpu-res/'%sector
#     fnames = os.listdir(res_dir)
#     for f in fnames:
#         cat = pd.read_csv(res_dir+f, header=None, skiprows=1)
#         cat_ticid = cat[0].to_numpy().astype('int')

#         ra= np.append(ra, cat[1])
#         dec = np.append(dec, cat[2])
#         pow = np.append(pow, cat[3])
#         snr = np.append(snr, cat[4])
#         dphi = np.append(dphi, cat[13])
#         wid = np.append(wid, np.int64(cat[5]))
#         nt = np.append(nt, np.int64(cat[12]))
#         per = np.append(per, cat[6])

#         suffix = f.split('.')[0][-4:]
#         plot_dir =  data_dir + 's%04d-ZTF/'%sector + 's%04d-bls'%sector+suffix+'/'
#         plot_fnames = np.array(os.listdir(plot_dir))
#         plot_ticid = np.array([int(plot.split('_')[8][3:]) for plot in plot_fnames])
#         for ticid in cat_ticid:
#             ind = np.nonzero(plot_ticid == ticid)[0][0]
#             fnames_jvr.append(plot_dir+plot_fnames[ind])

# fnames_jvr = np.array(fnames_jvr)


# # gmag = []
# # for i in range(len(ra)):
# #     coord = SkyCoord(ra=ra[i], dec=dec[i],
# #                      unit=(u.degree, u.degree), frame='icrs')
# #     try:
# #         width, height = u.Quantity(2, u.arcsec), u.Quantity(2, u.arcsec)
# #         r = Gaia.query_object_async(coordinate=coord, width=width, height=height)
# #         # j = Gaia.cone_search_async(co[i], radius=u.Quantity(3, u.arcsec))
# #     except:
# #         gmag.append(np.nan)
# #     if len(r['phot_g_mean_mag']) > 0:
# #         gmag.append(r['phot_g_mean_mag'][0])
# #     else:
# #         gmag.append(np.nan)
# # gmag = np.array(gmag)
# # np.savetxt(data_dir+'Gmag_JVR.txt', gmag)
# gmag = np.loadtxt(data_dir + 'Gmag_JVR.txt')

# fig4, ax4 = plt.subplots()
# ax4.plot(gmag_all, per_all[good_idx_all], '.k', ms=1,alpha=0.5)
# ax4.set_xlabel('Gaia Gmag')
# ax4.set_ylabel('Period [days]')
# ax4.plot(gmag[match], per[match], '<b', label='JVR-recovered')
# ax4.plot(gmag[~match], per_ztf[~match], 'Xr', label='JVR-unrecovered\nTrue period')
# ax4.plot(gmag[~match], per[~match], 'Xm', label='JVR-unrecovered\nBLS period')
# ax4.legend()
# fig4.savefig(out_dir+'gmag_per_JVR.png', dpi=300)
# ax4.set_ylim([0, 0.4])
# fig4.savefig(out_dir+'gmag_per_JVR_zoom.png', dpi=300)


# def get_stellar_density(ra, dec):
#     from astroquery.vizier import Vizier
#     import astropy.units as u

#     # Set the coordinates and search radius
#     coords = f"{ra} {dec}"
#     radius = 1 * u.arcmin # Search radius in arcminutes

#     # Query the VizieR database (using the Gaia catalog)
#     catalog = "I/345/gaia2"
#     v = Vizier(columns=['**'], catalog=catalog)
#     result = v.query_region(coords, radius=radius, catalog=catalog)

#     # Get the number of stars within the search radius
#     if len(result[0]) == 0:
#         star_count, density = np.nan, np.nan
#     else:
#         star_count = len(result[0])

#         # Calculate the stellar density (assuming uniform area)
#         density = star_count / (3.14 * radius**2)
#         density = density.value

#     return density

# # stellar_density = []
# # for i in range(len(ra)):
# #     stellar_density.append(get_stellar_density(ra[i], dec[i]))
# # stellar_density = np.array(stellar_density)
# # np.savetxt(data_dir+'stellar_density_JVR.txt', stellar_density)
# stellar_density = np.loadtxt(data_dir + 'stellar_density_JVR.txt')


# fig4, ax4 = plt.subplots()
# ax4.set_xlabel('Stellar density within 1 arcmin')
# ax4.set_ylabel('Period [days]')
# ax4.plot(stellar_density[match], per[match], '<b', label='JVR-recovered')
# ax4.plot(stellar_density[~match], per_ztf[~match], 'Xr', label='JVR-unrecovered\nTrue period')
# ax4.plot(stellar_density[~match], per[~match], 'Xm', label='JVR-unrecovered\nBLS period')
# ax4.legend()
# fig4.savefig(out_dir+'stellar_density_per_JVR.png', dpi=300)
# ax4.set_ylim([0, 0.4])
# fig4.savefig(out_dir+'stellar_density_per_JVR_zoom.png', dpi=300)

# # # >> inspect BLS plots for unrecovered JVR 
# # os.makedirs(out_dir+'unrecovered_JVR/', exist_ok=True)
# # for f in fnames_jvr[~match]:
# #     os.system('cp '+f+' '+out_dir+'unrecovered_JVR/')
# # plot_dirs = np.array([fn.split('/')[4] for fn in fnames_jvr[~match]])
# # labels, counts = np.unique(plot_dirs, return_counts=True)
# # os.makedirs(out_dir+'unrecovered_JVR_bls/', exist_ok=True)

# # result = []
# # for period, f in zip(per_ztf[~match], fnames_jvr[~match]):
# #     new_f = f.split('/')[-1][:-4] + '-true_period.png'
# #     ticid = np.int64(f.split('_')[8][3:])
# #     cam = np.int64(f.split('_')[11])
# #     ccd = np.int64(f.split('_')[13])
    

# #     ticid_ccd = np.load(lc_dir+'id-{}-{}.npy'.format(cam, ccd))
# #     ind = np.nonzero(ticid_ccd == ticid)[0]
# #     t = np.load(lc_dir+'ts-{}-{}.npy'.format(cam,ccd))
# #     y = np.load(lc_dir+'lc-{}-{}.npy'.format(cam,ccd))[ind]
# #     coord=np.load(lc_dir+'co-{}-{}.npy'.format(cam,ccd))[ind][0]
# #     cn = np.load(lc_dir+'cn-{}-{}.npy'.format(cam,ccd))

# #     t, y, cn = lcu.rm_qflag(t, y, cn, qflag_dir, sector, cam, ccd)
# #     t, y, flag = lcu.prep_lc(t, y)
# #     folded_t, folded_y, folded_dy = lcu.bin_timeseries(t%period, y, 100)
# #     fig, ax = plt.subplots()
# #     lcu.plot_phase_curve(ax, folded_t, folded_y, folded_dy, period=period)
# #     fig.savefig(out_dir+'unrecovered_JVR/'+new_f)
# #     print('Saved '+out_dir+'unrecovered_JVR/'+new_f)

# #     from Period_Finding import BLS
# #     dt = 5 # minutes
# #     pmin = period*1440. - dt # minutes
# #     pmax = period + dt/1440. # days
# #     dy = np.ones(y.shape)*0.1
# #     freqs_to_remove=lcu.rm_freq_tess()        
# #     t, y, dy, period, bls_power_best, freqs, power, q, phi0 = \
# #         BLS(t,y,dy,pmin=pmin,pmax=pmax,freqs_to_remove=freqs_to_remove)
# #     suffix = '_TIC%016d'%ticid+'_s%04d_'%sector+'cam_'+\
# #              str(cam)+'_ccd_'+str(ccd)+\
# #              '_ra_{}_dec_{}_'.format(coord[0], coord[1])
# #     res = lcu.vet_plot(t, y, freqs, power, q, phi0, output_dir=out_dir+'unrecovered_JVR_bls/',
# #                    objid=ticid, objid_type=None, suffix=suffix, 
# #                    ra=coord[0], dec=coord[1],
# #                    wd_main=wd_main, rp_ext=rp_ext, wd_tab=wd_tab)
# #     result.append([ticid, coord[0], coord[1]] + list(res))

# # np.savetxt(out_dir+'unrecovered_JVR_bls/GPU.result', np.array(result),
# #            fmt='%s,%10.5f,%10.5f,%10.5f,%10.5f,%i,%10.8f,%10.5f,%10.5f,%10.5f,%10.5f,%10.5f,%i,%10.5f',
# #            header='ticid, ra, dec, sig, snr, wid, period, period_min, q, phi0, epo, rp, nt, dphi')


# # plt.figure(figsize=(10,6))
# # x = range(len(labels))
# # plt.bar(x, counts)
# # _ = plt.xticks(x, labels, rotation='vertical')
# # plt.tight_layout() 
# # plt.savefig(out_dir+'unrecovered_JVR_bar.png')




# plt.figure()
# _=plt.hist(power_all[np.nonzero(power_all < 3000)], bins=100)
# plt.xlabel('Peak Significance')
# plt.ylabel('Num objects in TESS S61')
# plt.savefig(out_dir + 'pow_hist.png')
# print(out_dir + 'pow_hist.png')

# inds = np.nonzero( (power_all > 16) * (snr_all > 5) * (nt_all>130) )
# np.savetxt('/scratch/data/tess/lcur/ffi/cvae_data/metric_cut.txt', np.int64(ticid_all[inds]))

# # >> Inspect periodic sources taht aren't DWDs, UCBs or JVRs
# save_dir = out_dir + 'new_periodic/'
# os.makedirs(save_dir, exist_ok=True)
# co_wd = SkyCoord(ra=ra_all, dec=dec_all, unit=(u.degree, u.degree), frame='icrs')
# wd_idx = []

# co_dwd = SkyCoord(ra=ra_dwd, dec=dec_dwd, unit=(u.degree, u.degree), frame='icrs')
# idx, d2d, _ = co.match_to_catalog_sky(co_wd)
# good_idx = np.nonzero(d2d.to(u.arcsec).value > 2)[0] # want unmatched
# if len(good_idx) > 0:
#     wd_idx.extend(idx[good_idx])

# co_ucb = SkyCoord(ra=ra_ucb, dec=dec_ucb, unit=(u.degree, u.degree), frame='icrs')
# idx, d2d, _ = co.match_to_catalog_sky(co_wd)
# good_idx = np.nonzero(d2d.to(u.arcsec).value > 2)[0] # want unmatched
# if len(good_idx) > 0:
#     wd_idx.extend(idx[good_idx])

# co_jvr = SkyCoord(ra=ra, dec=dec, unit=(u.degree, u.degree), frame='icrs')
# idx, d2d, _ = co.match_to_catalog_sky(co_wd)
# good_idx = np.nonzero(d2d.to(u.arcsec).value > 2)[0] # want unmatched
# if len(good_idx) > 0:
#     wd_idx.extend(idx[good_idx])

# wd_idx = np.array(wd_idx)


# for cam in [1,2,3,4]:
#     for ccd in [1,2,3,4]:
#         suffix = '-{}-{}'.format(cam,ccd)
#         plot_dir =  data_dir + 's%04d-WD/'%sector + 's%04d-bls'%sector+suffix+'/'
#         plot_fnames = np.array(os.listdir(plot_dir))
#         plot_ticid = np.array([int(plot.split('_')[8][3:]) for plot in plot_fnames])
#         for f, ticid in zip(plot_fnames, plot_ticid):
#             if ticid in ticid_all[wd_idx]: # if unmatched
#                 os.system('cp '+plot_dir+f+' '+save_dir)
