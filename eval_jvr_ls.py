import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from astropy.coordinates import SkyCoord
import astropy.units as u 
from astroquery.gaia import Gaia
import lc_utils as lcu
Gaia.ROW_LIMIT = 5
Gaia.MAIN_GAIA_TABLE = "gaiadr3.gaia_source"

sector_dir = "/scratch/echickle/s0056_LS/"
stat_dir = sector_dir + "stats/"
vet_dir = sector_dir + "vet/"
os.makedirs(vet_dir, exist_ok=True)

lc_dir = '/scratch/data/tess/lcur/ffi/s0056-lc-ZTF/'
qflag_dir = "/scratch/echickle/QLPqflags/"
wd_tab= "/scratch/echickle/WDs.txt"
wd_main = "/scratch/echickle/GaiaEDR3_WD_main.fits"
rp_ext = "/scratch/echickle/GaiaEDR3_WD_RPM_ext.fits"
sector = 56

pmax = 0.15*1440
cat = np.loadtxt("/scratch/echickle/ZTF_Eclipses.txt", usecols=(1,2,3))
ra_ztf, dec_ztf = cat[:,0], cat[:,1]
per_ztf = cat[:,2]

wd_cat  = pd.read_csv(wd_tab, header=None, sep='\s+', dtype='str')
ra_wd = wd_cat[1].to_numpy().astype('float')
dec_wd = wd_cat[2].to_numpy().astype('float')
gmag_wd = wd_cat[4].to_numpy().astype('float')

# # -- TESS signals --------------------------------------------------------------

# # >> all signals
# ticid_all = []
# ra_all = []
# dec_all = []
# power_all = []
# snr_all = []
# wid_all = []
# per_all = []
# nt_all = []
# dphi_all = []

# # >> signals discovered in S61
# ticid_tess = np.array([803489769, 36085812, 800042858, 270679102, 455206965, 452954413, 767706310, 96852226, 677736827])
# ra_tess = []
# dec_tess = []
# pow_tess = []
# snr_tess = []
# wid_tess = [d]
# per_tess = []
# nt_tess = []
# dphi_tess= []

# # >> KB UCBs systems in the white dwarf catalog
# ticid_ucb = np.array([719830487, 702774311, 713144139, 1939638541, 746047312, 1504732386, 1717410027, 1958881490, 2040677137, 1270451609])
# per_ucb_true = np.array([0.010030088, 0.014240466, 0.014282094, 0.0144914, 0.016464692, 0.018356874, 0.0281957, 0.029969971, 0.038365738, 0.04542631929])
# ra_ucb = []
# dec_ucb = []
# pow_ucb = []
# snr_ucb = []
# wid_ucb = []
# per_ucb = []
# nt_ucb = []
# dphi_ucb = []
# match_ucb = []
# per_ucb_match = []
# tic_ucb_match = []
# no_match_ucb = []

# # >> long period DWDs in the white dwarf catalog
# ticid_dwd = np.array([755049150, 903242599, 840220096, 219868627])
# per_dwd_true = np.array([0.0800126, 0.09986, 0.11601549, 0.246137])
# ra_dwd = []
# dec_dwd = []
# pow_dwd = []
# snr_dwd = []
# wid_dwd = []
# per_dwd = []
# nt_dwd = []
# dphi_dwd = []
# match_dwd = []
# per_dwd_match = []
# tic_dwd_match = []
# no_match_dwd = []

# # ticid, ra, dec, sig, snr, wid, period, period_min, q, phi0, dur, epo, rp, nt, dphi
# # # ticid, ra, dec, sig, snr, wid, period, period_min, q, phi0, epo, rp, nt, dphi
# for sector in [56]: # !!
#     for cam in [1,2,3,4]:
#         for ccd in [1,2,3,4]: 
#             f = '/scratch/echickle/s%04d-WD/'%sector+'s%04d-gpu-res/'%sector+\
#             'GPU-{}-{}-{}.result'.format(sector, cam, ccd)
#             cat = pd.read_csv(f, header=None, skiprows=1)
#             # ra_all.extend(cat[1])
#             # dec_all.extend(cat[2])
#             # power_all.extend(cat[3])
#             # snr_all.extend(cat[4])
#             # wid_all.extend(np.int64(cat[5]))
#             # per_all.extend(cat[6])
#             # nt_all.extend(np.int64(cat[13]))
#             # dphi_all.extend(cat[14])
#             ticid_all.extend(cat[0])
#             ra_all.extend(cat[1])
#             dec_all.extend(cat[2])
#             power_all.extend(cat[3])
#             snr_all.extend(cat[4])
#             wid_all.extend(np.int64(cat[5]))
#             per_all.extend(cat[6])
#             nt_all.extend(np.int64(cat[12]))
#             dphi_all.extend(cat[13])


#             cat_ticid = np.int64(cat[0])
#             _, inds, inds1 = np.intersect1d(cat_ticid, ticid_tess, return_indices=True)
#             ra_tess.extend(cat[1][inds])
#             dec_tess.extend(cat[2][inds])
#             pow_tess.extend(cat[3][inds])
#             snr_tess.extend(cat[4][inds])
#             wid_tess.extend(np.int64(cat[5])[inds])
#             per_tess.extend(cat[6][inds])
#             nt_tess.extend(np.int64(cat[12])[inds])
#             dphi_tess.extend(cat[13][inds])

#             _, inds, inds1 = np.intersect1d(cat_ticid, ticid_ucb, return_indices=True)
#             ra_ucb.extend(cat[1][inds])
#             dec_ucb.extend(cat[2][inds])
#             pow_ucb.extend(cat[3][inds])
#             snr_ucb.extend(cat[4][inds])
#             wid_ucb.extend(np.int64(cat[5])[inds])
#             per_bls = cat[6][inds].to_numpy()
#             per_ucb.extend(per_bls)
#             nt_ucb.extend(np.int64(cat[12])[inds])
#             dphi_ucb.extend(cat[13][inds])

#             dt = 2/1440
#             for i in range(len(per_bls)):
#                 if (per_bls[i] > per_ucb_true[inds1][i]-dt \
#                     and per_bls[i] < per_ucb_true[inds1][i]+dt) or\
#                 (per_bls[i] > per_ucb_true[inds1][i]*2-dt \
#                  and per_bls[i] < per_ucb_true[inds1][i]*2+dt) or\
#                 (per_bls[i] > per_ucb_true[inds1][i]/2-dt \
#                  and per_bls[i] < per_ucb_true[inds1][i]/2+dt):
#                     match_ucb.append(True)
#                 else:
#                     match_ucb.append(False)
#                 per_ucb_match.append(per_ucb_true[inds1][i])
#                 tic_ucb_match.append(ticid_ucb[inds1][i])

#             _, inds, inds1 = np.intersect1d(cat_ticid, ticid_dwd, return_indices=True)
#             ra_dwd.extend(cat[1][inds])
#             dec_dwd.extend(cat[2][inds])
#             pow_dwd.extend(cat[3][inds])
#             snr_dwd.extend(cat[4][inds])
#             wid_dwd.extend(np.int64(cat[5])[inds])
#             per_bls = cat[6][inds].to_numpy()
#             per_dwd.extend(per_bls)
#             nt_dwd.extend(np.int64(cat[12])[inds])
#             dphi_dwd.extend(cat[13][inds])

#             dt = 2/1440
#             for i in range(len(per_bls)):
#                 if (per_bls[i] > per_dwd_true[inds1][i]-dt \
#                     and per_bls[i] < per_dwd_true[inds1][i]+dt) or\
#                 (per_bls[i] > per_dwd_true[inds1][i]*2-dt \
#                  and per_bls[i] < per_dwd_true[inds1][i]*2+dt) or\
#                 (per_bls[i] > per_dwd_true[inds1][i]/2-dt \
#                  and per_bls[i] < per_dwd_true[inds1][i]/2+dt):
#                     match_dwd.append(True)
#                 else:
#                     # no_match_dwd.append(inds1[i])
#                     match_dwd.append(False)
#                 per_dwd_match.append(per_dwd_true[inds1][i])
#                 tic_dwd_match.append(ticid_dwd[inds1][i])

# ticid_all = np.array(ticid_all)
# power_all = np.array(power_all)
# snr_all = np.array(snr_all)
# wid_all = np.array(wid_all)
# per_all = np.array(per_all)
# nt_all = np.array(nt_all)
# dphi_all = np.array(dphi_all)

# pow_ucb = np.array(pow_ucb)
# snr_ucb = np.array(snr_ucb)
# wid_ucb = np.array(wid_ucb)
# per_ucb = np.array(per_ucb)
# nt_ucb = np.array(nt_ucb)
# dphi_ucb = np.array(dphi_ucb)
# match_ucb = np.array(match_ucb)
# per_ucb_match = np.array(per_ucb_match)
# tic_ucb_match = np.array(tic_ucb_match)

# pow_dwd = np.array(pow_dwd)
# snr_dwd = np.array(snr_dwd)
# wid_dwd = np.array(wid_dwd)
# per_dwd = np.array(per_dwd)
# nt_dwd = np.array(nt_dwd)
# dphi_dwd = np.array(dphi_dwd)
# match_dwd = np.array(match_dwd)
# per_dwd_match = np.array(per_dwd_match)
# tic_dwd_match = np.array(tic_dwd_match)

# -- initialize metric plots ---------------------------------------------------

fig1, ax1 = plt.subplots()
# ax1.plot(power_all, snr_all, '.k', ms=1, alpha=0.5)
# ax1.plot(pow_tess, snr_tess, 'vr', label='TESS')

fig2, ax2 = plt.subplots()
# ax2.plot(power_all, wid_all, '.k', ms=1, alpha=0.5)
# ax2.plot(pow_tess, wid_tess, 'vr', label='TESS')

fig3, ax3 = plt.subplots()
# ax3.plot(nt_all, dphi_all, '.k', ms=1, alpha=0.5)
# ax3.plot(nt_tess, dphi_tess, 'vr', label='TESS')

# if len(match_ucb) > 0:
#     ax1.plot(pow_ucb[match_ucb], snr_ucb[match_ucb], '>g', label='KBUCB')
#     ax2.plot(pow_ucb[match_ucb], wid_ucb[match_ucb], '>g', label='KBUCB')
#     ax3.plot(nt_ucb[match_ucb], dphi_ucb[match_ucb], '>g', label='KBUCB')
# if len(match_dwd) > 0:
#     ax1.plot(pow_dwd[match_dwd], snr_dwd[match_dwd], 'oc', label='DWD')
#     ax2.plot(pow_dwd[match_dwd], wid_dwd[match_dwd], 'oc', label='DWD')
#     ax3.plot(nt_dwd[match_dwd], dphi_dwd[match_dwd], 'oc', label='DWD')

# -- recovery plots ------------------------------------------------------------


# print('KB UCB Recovery: '+str(np.count_nonzero(match_ucb))+' / '+str(len(per_ucb)))
# print('Recovered UCBs: '+','.join(tic_ucb_match.astype('str')))
# print('DWD Recovery: '+str(np.count_nonzero(match_dwd))+' / '+str(len(per_dwd)))
# print('Recovered DWDs: '+','.join(tic_dwd_match.astype('str')))

# co_wd = SkyCoord(ra=ra_wd, dec=dec_wd, unit=(u.degree, u.degree), frame='icrs')

# co = SkyCoord(ra=ra_all, dec=dec_all, unit=(u.degree, u.degree), frame='icrs')
# idx, d2d, _ = co.match_to_catalog_sky(co_wd)
# good_idx_all = np.nonzero(d2d.to(u.arcsec).value < 2)
# gmag_all = gmag_wd[idx[good_idx_all]]

# co = SkyCoord(ra=ra_ucb, dec=dec_ucb, unit=(u.degree, u.degree), frame='icrs')
# idx, d2d, _ = co.match_to_catalog_sky(co_wd)
# good_idx = np.nonzero(d2d.to(u.arcsec).value < 2)
# gmag = gmag_wd[idx[good_idx]]

# fig4, ax4 = plt.subplots()
# ax4.plot(gmag_all, per_all[good_idx_all], '.k', ms=1,alpha=0.5)
# ax4.set_xlabel('Gaia Gmag')
# ax4.set_ylabel('Period [days]')
# if len(match_ucb) > 0:
#     ax4.plot(gmag[match_ucb], per_ucb[match_ucb], '>g', label='KBUCB-recovered')
#     ax4.plot(gmag[~match_ucb], per_ucb_match[~match_ucb], 'Xr', label='KBUCB-unrecovered\nTrue period')
#     ax4.plot(gmag[~match_ucb], per_ucb[~match_ucb], 'Xm', label='KBUCB-unrecovered\nBLS period')
# ax4.legend()
# fig4.savefig(out_dir+'gmag_per_KBUCB.png', dpi=300)
# ax4.set_ylim([0, 0.4])
# fig4.savefig(out_dir+'gmag_per_KBUCB_zoom.png', dpi=300)

# co = SkyCoord(ra=ra_dwd, dec=dec_dwd, unit=(u.degree, u.degree), frame='icrs')
# co_wd = SkyCoord(ra=ra_wd, dec=dec_wd, unit=(u.degree, u.degree), frame='icrs')
# idx, d2d, _ = co.match_to_catalog_sky(co_wd)
# good_idx = idx[np.nonzero(d2d.to(u.arcsec).value < 2)]
# gmag = gmag_wd[good_idx]

# fig4, ax4 = plt.subplots()
# ax4.plot(gmag_all, per_all[good_idx_all], '.k', ms=1,alpha=0.5)
# ax4.set_xlabel('Gaia Gmag')
# ax4.set_ylabel('Period [days]')
# if len(match_dwd) > 0:
#     ax4.plot(gmag[match_dwd], per_dwd[match_dwd], 'oc', label='Long DWD-recovered')
#     ax4.plot(gmag[~match_dwd], per_dwd_match[~match_dwd], 'Xr', label='Long DWD-unrecovered\nTrue period')
#     ax4.plot(gmag[~match_dwd], per_dwd[~match_dwd], 'Xm', label='Long DWD-unrecovered\nBLS period')
# ax4.legend()
# fig4.savefig(out_dir+'gmag_per_DWD.png', dpi=300)
# ax4.set_ylim([0, 0.4])
# fig4.savefig(out_dir+'gmag_per_DWD_zoom.png', dpi=300)


# -- JVR -----------------------------------------------------------------------

# ticid
ra, dec, sig, snr, dphi, wid, per = [np.empty(0)]*7
fnames_jvr = []

# ticid, ra, dec, sig, wid, per, per_min, dphi
for sector in [56]: # !!
    fnames = os.listdir(stat_dir)
    for f in fnames:
        cat = pd.read_csv(stat_dir+f, header=None, skiprows=1)
        cat_ticid = cat[0].to_numpy().astype('int')

        ra= np.append(ra, cat[1])
        dec = np.append(dec, cat[2])
        sig = np.append(sig, cat[3])
        wid = np.append(wid, np.int64(cat[4]))
        per = np.append(per, cat[5])
        dphi = np.append(dphi, cat[7])

        cam = int(f.split('-')[2])
        ccd = int(f.split('-')[3][:1])
        plot_dir =  sector_dir + 'cam{}-ccd{}/'.format(cam,ccd)
        plot_fnames = np.array(os.listdir(plot_dir))
        plot_ticid = np.array([int(plot.split('_')[6][3:]) for plot in plot_fnames])
        for ticid in cat_ticid:
            ind = np.nonzero(plot_ticid == ticid)[0][0]
            fnames_jvr.append(plot_dir+plot_fnames[ind])

fnames_jvr = np.array(fnames_jvr)

co = SkyCoord(ra=ra, dec=dec, unit=(u.degree, u.degree), frame='icrs')
co_ztf = SkyCoord(ra=ra_ztf, dec=dec_ztf, unit=(u.degree, u.degree), frame='icrs')
idx, d2d, _ = co.match_to_catalog_sky(co_ztf)
good_idx = idx[np.nonzero(d2d.to(u.arcsec).value < 2)]
per_ztf = per_ztf[good_idx]

match = []
no_match = []
dt = 2/1440
for i in range(len(per)):
    if (per[i] > per_ztf[i]-dt and per[i] < per_ztf[i]+dt) or\
    (per[i] > per_ztf[i]*2-dt and per[i] < per_ztf[i]*2+dt) or\
    (per[i] > per_ztf[i]/2-dt and per[i] < per_ztf[i]/2+dt):
        match.append(True)
    else:
        # no_match.append(i)
        match.append(False)
match, no_match = np.array(match), np.array(no_match)

print('JVR Recovery: '+str(np.count_nonzero(match))+' / '+str(len(per)))

# res = np.array([ra_list, dec_list, loc, power, snr, wid, per, per_ztf_list, match], dtype='str').T     
# np.savetxt(out_dir+'eval_jvr.txt', res, fmt='%s',
#            header='RA Dec Pow SNR Wid Per Per_ZTF Match')

# ax1.plot(sig[match], snr[match], '<b', label='JVR')
ax2.plot(sig[match], wid[match], '<b', label='JVR')
# ax3.plot(nt[match], dphi[match], '<b', label='JVR')

# gmag = []
# for i in range(len(ra)):
#     coord = SkyCoord(ra=ra[i], dec=dec[i],
#                      unit=(u.degree, u.degree), frame='icrs')
#     try:
#         width, height = u.Quantity(2, u.arcsec), u.Quantity(2, u.arcsec)
#         r = Gaia.query_object_async(coordinate=coord, width=width, height=height)
#         # j = Gaia.cone_search_async(co[i], radius=u.Quantity(3, u.arcsec))
#     except:
#         gmag.append(np.nan)
#     if len(r['phot_g_mean_mag']) > 0:
#         gmag.append(r['phot_g_mean_mag'][0])
#     else:
#         gmag.append(np.nan)
# gmag = np.array(gmag)
# np.savetxt(vet_dir+'Gmag_JVR.txt', gmag)
gmag = np.loadtxt(vet_dir + 'Gmag_JVR.txt')

fig4, ax4 = plt.subplots()
# ax4.plot(gmag_all, per_all[good_idx_all], '.k', ms=1,alpha=0.5)
ax4.set_xlabel('Gaia Gmag')
ax4.set_ylabel('Period [days]')
ax4.plot(gmag[match], per[match], '<b', label='JVR-recovered')
ax4.plot(gmag[~match], per_ztf[~match], 'Xr', label='JVR-unrecovered\nTrue period')
ax4.plot(gmag[~match], per[~match], 'Xm', label='JVR-unrecovered\nBLS period')
ax4.legend()
fig4.savefig(vet_dir+'gmag_per_JVR.png', dpi=300)
ax4.set_ylim([0, 0.4])
fig4.savefig(vet_dir+'gmag_per_JVR_zoom.png', dpi=300)


def get_stellar_density(ra, dec):
    from astroquery.vizier import Vizier
    import astropy.units as u

    # Set the coordinates and search radius
    coords = f"{ra} {dec}"
    radius = 1 * u.arcmin # Search radius in arcminutes

    # Query the VizieR database (using the Gaia catalog)
    catalog = "I/345/gaia2"
    v = Vizier(columns=['**'], catalog=catalog)
    result = v.query_region(coords, radius=radius, catalog=catalog)

    # Get the number of stars within the search radius
    if len(result[0]) == 0:
        star_count, density = np.nan, np.nan
    else:
        star_count = len(result[0])

        # Calculate the stellar density (assuming uniform area)
        density = star_count / (3.14 * radius**2)
        density = density.value

    return density

# stellar_density = []
# for i in range(len(ra)):
#     stellar_density.append(get_stellar_density(ra[i], dec[i]))
# stellar_density = np.array(stellar_density)
# np.savetxt(vet_dir+'stellar_density_JVR.txt', stellar_density)
stellar_density = np.loadtxt(vet_dir + 'stellar_density_JVR.txt')


fig4, ax4 = plt.subplots()
ax4.set_xlabel('Stellar density within 1 arcmin')
ax4.set_ylabel('Period [days]')
ax4.plot(stellar_density[match], per[match], '<b', label='JVR-recovered')
ax4.plot(stellar_density[~match], per_ztf[~match], 'Xr', label='JVR-unrecovered\nTrue period')
ax4.plot(stellar_density[~match], per[~match], 'Xm', label='JVR-unrecovered\nBLS period')
ax4.legend()
fig4.savefig(vet_dir+'stellar_density_per_JVR.png', dpi=300)
ax4.set_ylim([0, 0.4])
fig4.savefig(vet_dir+'stellar_density_per_JVR_zoom.png', dpi=300)

# # >> inspect BLS plots for unrecovered JVR 
# os.makedirs(out_dir+'unrecovered_JVR/', exist_ok=True)
# for f in fnames_jvr[~match]:
#     os.system('cp '+f+' '+out_dir+'unrecovered_JVR/')
# plot_dirs = np.array([fn.split('/')[4] for fn in fnames_jvr[~match]])
# labels, counts = np.unique(plot_dirs, return_counts=True)
# os.makedirs(out_dir+'unrecovered_JVR_bls/', exist_ok=True)

# result = []
# for period, f in zip(per_ztf[~match], fnames_jvr[~match]):
#     new_f = f.split('/')[-1][:-4] + '-true_period.png'
#     ticid = np.int64(f.split('_')[8][3:])
#     cam = np.int64(f.split('_')[11])
#     ccd = np.int64(f.split('_')[13])
    

#     ticid_ccd = np.load(lc_dir+'id-{}-{}.npy'.format(cam, ccd))
#     ind = np.nonzero(ticid_ccd == ticid)[0]
#     t = np.load(lc_dir+'ts-{}-{}.npy'.format(cam,ccd))
#     y = np.load(lc_dir+'lc-{}-{}.npy'.format(cam,ccd))[ind]
#     coord=np.load(lc_dir+'co-{}-{}.npy'.format(cam,ccd))[ind][0]
#     cn = np.load(lc_dir+'cn-{}-{}.npy'.format(cam,ccd))

#     t, y, cn = lcu.rm_qflag(t, y, cn, qflag_dir, sector, cam, ccd)
#     t, y, flag = lcu.prep_lc(t, y)
#     folded_t, folded_y, folded_dy = lcu.bin_timeseries(t%period, y, 100)
#     fig, ax = plt.subplots()
#     lcu.plot_phase_curve(ax, folded_t, folded_y, folded_dy, period=period)
#     fig.savefig(out_dir+'unrecovered_JVR/'+new_f)
#     print('Saved '+out_dir+'unrecovered_JVR/'+new_f)

#     from Period_Finding import BLS
#     dt = 5 # minutes
#     pmin = period*1440. - dt # minutes
#     pmax = period + dt/1440. # days
#     dy = np.ones(y.shape)*0.1
#     freqs_to_remove=lcu.rm_freq_tess()        
#     t, y, dy, period, bls_power_best, freqs, power, q, phi0 = \
#         BLS(t,y,dy,pmin=pmin,pmax=pmax,freqs_to_remove=freqs_to_remove)
#     suffix = '_TIC%016d'%ticid+'_s%04d_'%sector+'cam_'+\
#              str(cam)+'_ccd_'+str(ccd)+\
#              '_ra_{}_dec_{}_'.format(coord[0], coord[1])
#     res = lcu.vet_plot(t, y, freqs, power, q, phi0, output_dir=out_dir+'unrecovered_JVR_bls/',
#                    objid=ticid, objid_type=None, suffix=suffix, 
#                    ra=coord[0], dec=coord[1],
#                    wd_main=wd_main, rp_ext=rp_ext, wd_tab=wd_tab)
#     result.append([ticid, coord[0], coord[1]] + list(res))

# np.savetxt(out_dir+'unrecovered_JVR_bls/GPU.result', np.array(result),
#            fmt='%s,%10.5f,%10.5f,%10.5f,%10.5f,%i,%10.8f,%10.5f,%10.5f,%10.5f,%10.5f,%10.5f,%i,%10.5f',
#            header='ticid, ra, dec, sig, snr, wid, period, period_min, q, phi0, epo, rp, nt, dphi')


# plt.figure(figsize=(10,6))
# x = range(len(labels))
# plt.bar(x, counts)
# _ = plt.xticks(x, labels, rotation='vertical')
# plt.tight_layout() 
# plt.savefig(out_dir+'unrecovered_JVR_bar.png')

# -- save figures --------------------------------------------------------------

ax1.set_xlabel('Peak Significance = (peak - median)/MAD')
ax1.set_ylabel('SNR = depth/MAD')
ax1.legend()
fig1.tight_layout()
fig1.savefig(vet_dir + 'pow_snr.png')
print(vet_dir + 'pow_snr.png')

# ax1.set_xlim([0, 200])
# ax1.set_ylim([0, max(snr_tess)])
# fig1.savefig(out_dir + 'pow_snr_zoom.png')
# print(out_dir + 'pow_snr_zoom.png')

ax2.set_xlabel('Peak Significance = (peak - median)/MAD')
ax2.set_ylabel('Peak width')
ax2.legend()
fig2.tight_layout()
fig2.savefig(vet_dir + 'pow_wid.png')
print(vet_dir + 'pow_wid.png')

# ax2.set_xlim([0, 200])
# ax2.set_ylim([0, max(wid_tess)])
# fig2.savefig(out_dir + 'pow_wid_zoom.png')
# print(out_dir + 'pow_wid_zoom.png')

ax3.set_xlabel('Number of points in transit')
ax3.set_ylabel('Dphi')
ax3.legend()
fig3.tight_layout()
fig3.savefig(vet_dir + 'nt_dphi.png')
print(vet_dir + 'nt_dphi.png')

ax3.set_ylim([0, 0.04])
fig3.savefig(vet_dir + 'nt_dphi_zoom.png')

os.makedirs(vet_dir+'recovered_JVR/', exist_ok=True)
for i in range(len(fnames_jvr[match])):
    os.system('cp '+fnames_jvr[match][i]+' '+vet_dir+'recovered_JVR/')

# plt.figure()
# _=plt.hist(power_all[np.nonzero(power_all < 3000)], bins=100)
# plt.xlabel('Peak Significance')
# plt.ylabel('Num objects in TESS S61')
# plt.savefig(out_dir + 'pow_hist.png')
# print(out_dir + 'pow_hist.png')

# inds = np.nonzero( (power_all > 16) * (snr_all > 5) * (nt_all>130) )
# np.savetxt('/scratch/data/tess/lcur/ffi/cvae_data/metric_cut.txt', np.int64(ticid_all[inds]))

# >> Inspect periodic sources taht aren't DWDs, UCBs or JVRs
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
