import numpy as np
import os

def read_result_file(fname, bls=True):
    import pandas as pd
    if bls:
        # ticid, ra, dec, sig, snr, wid, period, period_min, q, phi0, dur, epo, rp, nt, dphi
        cat = pd.read_csv(fname, header=None, skiprows=1)    
        ticid, ra, dec, power, snr, wid, per, nt, dphi = cat[0], cat[1], cat[2], cat[3], cat[4], cat[5], cat[6], cat[12], cat[13]
        return ticid, ra, dec, power, snr, wid, per, nt, dphi
    else:
        # ticid, ra, dec, sig, wid, per, per_min, dphi
        cat = pd.read_csv(fname, header=None, skiprows=1)    
        ticid, ra, dec, power, wid, per, dphi = cat[0], cat[1], cat[2], cat[3], cat[4], cat[5], cat[7]
        return ticid, ra, dec, power, wid, per, dphi
        

def append_result_file(result_dir, sector_list, cam_list=[1,2,3,4], ccd_list=[1,2,3,4], bls=True):
    if bls:
        result_list = [np.empty(0)]*9
    else:
        result_list = [np.empty(0)]*7
    for sector in sector_list:
        for cam in cam_list:
            for ccd in ccd_list:
                if bls:
                    fname = result_dir+'BLS-{}-{}-{}.result'.format(sector,cam,ccd)
                else:
                    fname = result_dir+'LS-{}-{}-{}.result'.format(sector,cam,ccd)
                if os.path.exists(fname):
                    result = read_result_file(fname, bls=bls)

                    # if np.count_nonzero( np.int64(result[0].to_numpy() ) == 0 ) >0:
                    #     import pdb
                    #     pdb.set_trace()

                    for i in range(len(result)):
                        result_list[i] = np.append(result_list[i], result[i])
    result_list = np.array(result_list)
    result_list = result_list.T
    return result_list
    
def match_period(period, period_true, dt=2./1440):
    if np.isscalar(period):
        period, period_true = [period], [period_true]
    
    match_list = []
    for i in range(len(period)):
        
        if (period[i] > period_true[i]-dt \
            and period[i] < period_true[i]+dt) or\
        (period[i] > period_true[i]*2-dt \
         and period[i] < period_true[i]*2+dt) or\
        (period[i] > period_true[i]/2-dt \
         and period[i] < period_true[i]/2+dt):
            match_list.append(True)
        else:
            match_list.append(False)
            
    if len(period) == 1:
        return match_list[0]
    else:
        return np.array(match_list)
    
def match_catalog(result_list, ticid_catalog, period_true=None):
    ticid = np.int64(result_list[:,0])
    _, inds, inds1 = np.intersect1d(ticid, ticid_catalog, return_indices=True)
    result_catalog = [result[inds] for result in result_list.T]
    result_catalog = np.array(result_catalog).T
    ticid_catalog = ticid_catalog[inds1]

    if period_true is not None:
        period_true = period_true[inds1]
        return ticid_catalog, period_true, result_catalog
    
    else:
        return ticid_catalog, result_catalog

def match_coord(ra_catalog, dec_catalog, per_catalog_true, result_list):
# co_dwd = SkyCoord(ra=ra_dwd, dec=dec_dwd, unit=(u.degree, u.degree), frame='icrs')
# idx, d2d, _ = co.match_to_catalog_sky(co_wd)
# good_idx = np.nonzero(d2d.to(u.arcsec).value > 2)[0] # want unmatched
# if len(good_idx) > 0:
#     wd_idx.extend(idx[good_idx])
    import astropy.units as u
    from astropy.coordinates import SkyCoord
    co_result = SkyCoord(ra=result_list[:,1], dec=result_list[:,2], unit=(u.degree, u.degree), frame='icrs')
    co_catalog = SkyCoord(ra=ra_catalog, dec=dec_catalog, unit=(u.degree, u.degree), frame='icrs')
    idx, d2d, _ = co_catalog.match_to_catalog_sky(co_result)
    good_idx = np.nonzero(d2d.to(u.arcsec).value < 2)[0]

    ra_catalog = ra_catalog[good_idx]
    dec_catalog = dec_catalog[good_idx]
    per_catalog_true = per_catalog_true[good_idx]
    result_list = result_list[idx[good_idx]]
    
    return ra_catalog, dec_catalog, per_catalog_true, result_list
    

def plot_gmag(out_dir, wd_tab, result_list, result_catalog, match_catalog, per_catalog_true, suffix='', gmag_catalog=None, y_max=2.):
    import matplotlib.pyplot as plt
    import pandas as pd
    import astropy.units as u
    from astropy.coordinates import SkyCoord
    
    fig, ax = plt.subplots(figsize=(5,4))
    ax.set_xlabel('Gaia Gmag')
    ax.set_ylabel('Period [days]')    
    
    # Load Gaia white darf catalog
    wd_cat  = pd.read_csv(wd_tab, header=None, sep='\s+', dtype='str')
    ra_wd = wd_cat[1].to_numpy().astype('float')
    dec_wd = wd_cat[2].to_numpy().astype('float')
    gmag_wd = wd_cat[4].to_numpy().astype('float')
    co_wd = SkyCoord(ra=ra_wd, dec=dec_wd, unit=(u.degree, u.degree), frame='icrs')
    
    # Match to result list (all results)
    co = SkyCoord(ra=result_list[:,1], dec=result_list[:,2], unit=(u.degree, u.degree), frame='icrs')
    idx, d2d, _ = co.match_to_catalog_sky(co_wd)
    good_idx_all = np.nonzero(d2d.to(u.arcsec).value < 2)
    gmag_all = gmag_wd[idx[good_idx_all]]
    per_all = result_list[:,6]
    ax.plot(gmag_all, per_all[good_idx_all], '.k', ms=1, alpha=0.5)
    
    if gmag_catalog is None:
        # Match to catalog
        co = SkyCoord(ra=result_catalog[:,1], dec=result_catalog[:,2], unit=(u.degree, u.degree), frame='icrs')
        idx, d2d, _ = co.match_to_catalog_sky(co_wd)
        good_idx = np.nonzero(d2d.to(u.arcsec).value < 2)
        gmag_catalog = gmag_wd[idx[good_idx]]
    
    if np.count_nonzero(match_catalog) > 0:
        per_catalog = result_catalog[:,6]
        ax.plot(gmag_catalog[~match_catalog], per_catalog_true[~match_catalog], 'Xr', label=suffix+'-unrecovered\nTrue period')
        ax.plot(gmag_catalog[~match_catalog], per_catalog[~match_catalog], 'Xm', label=suffix+'-unrecovered\nBLS period')
        ax.plot(gmag_catalog[match_catalog], per_catalog[match_catalog], '>g', label=suffix+'-recovered')        
    ax.legend()
    fig.tight_layout()
    fig.savefig(out_dir+'gmag_per_'+suffix+'.png', dpi=300)
    ax.set_ylim([0, y_max])
    fig.tight_layout()
    fig.savefig(out_dir+'gmag_per_'+suffix+'_zoom.png', dpi=300)
        
def get_gmag(result_list, mydir, overwrite=False):
    import astropy.units as u
    from astropy.coordinates import SkyCoord  
    import os
    from astroquery.gaia import Gaia
    
    fname = mydir+'Gmag_JVR.txt'
    if os.path.exists(fname) and not overwrite:
        gmag = np.loadtxt(fname)
    else:
        gmag = []
        ra, dec = result_list[:,1], result_list[:,2]
        for i in range(len(ra)):
            coord = SkyCoord(ra=ra[i], dec=dec[i],
                              unit=(u.degree, u.degree), frame='icrs')
            try:
                width, height = u.Quantity(2, u.arcsec), u.Quantity(2, u.arcsec)
                r = Gaia.query_object_async(coordinate=coord, width=width, height=height)
                # j = Gaia.cone_search_async(co[i], radius=u.Quantity(3, u.arcsec))
            except:
                gmag.append(np.nan)
            if len(r['phot_g_mean_mag']) > 0:
                gmag.append(r['phot_g_mean_mag'][0])
            else:
                gmag.append(np.nan)
        gmag = np.array(gmag)
        np.savetxt(fname, gmag)
    return gmag
