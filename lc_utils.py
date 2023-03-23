import numpy as np
import matplotlib.pyplot as plt
import pdb

def extract_lc():
    import pandas as pd
    
    # >> load white dwarf catalog
    wd_cat  = pd.read_csv('/home/echickle/work/tess_dwd/WDs.txt', header=None, sep='\s+')

def atlas_lc(f):
    data=np.loadtxt(f,usecols=(0,3,4,16),skiprows=0)    
    # data=np.loadtxt(f,usecols=(0,3,4,16),skiprows=0)
    # coords=np.loadtxt(f,usecols=(8,9),skiprows=0)
    # RA=np.mean(coords[:,0])
    # Dec=np.mean(coords[:,1])

    # >> filter by limiting magnitude
    fil=np.loadtxt(f,usecols=(5),skiprows=0,dtype=str)
    fil=fil[data[:,3]>17.5]
    data=data[data[:,3]>17.5]

    # return data, RA, Dec
    return data

def load_hipercam(logfile):
    '''Loads logfile produced by HiPERCAM reduce pipeline. 
    5 CCDs: u, g, r, i, z bands.'''

    # >> get column names
    with open(logfile, 'r') as f:
        cols = np.array(f.readlines()[417][6:-2].split(' '))

    data = np.loadtxt(logfile, skiprows=439)
    t = data[:,np.nonzero(cols == 'MJD')[0][0]]
    u = data[:,np.nonzero(cols == 'counts_1')[0][0]] # no
    g = data[:,np.nonzero(cols == 'counts_2')[0][0]]
    r = data[:,np.nonzero(cols == 'counts_3')[0][0]]

    return

def bin_timeseries(t, y, bins, dy=None):
    trunc = len(y) // bins * bins
    t_binned = np.split(t[:trunc], bins)
    y_binned = np.split(y[:trunc], bins)
    if type(dy) != type(None):
        dy_binned = np.split(dy[:trunc], bins)

    if trunc < len(y):
        y_binned[-1] = np.append(y_binned[-1], y[trunc:])
        t_binned[-1] = np.append(t_binned[-1], t[trunc:])
        if type(dy) != type(None):
            dy_binned[-1] = np.append(dy_binned[-1], dy[trunc:])

    if type(dy) == type(None):
        dy_binned = []
    for i in range(bins):
        t_binned[i] = np.average(t_binned[i])

        # dy_binned[i] = 1 / np.sqrt(np.sum(dy_binned[i]))
        npts = len(y_binned[i])
        
        if type(dy) == type(None):
            err = np.std(y_binned[i]) / np.sqrt(npts)            
            y_binned[i] = np.average(y_binned[i])
            dy_binned.append(err)
        else:
            y_binned[i] = np.average(y_binned[i], weights=1./dy_binned[i]**2)
            dy_binned[i] = np.average(dy_binned[i]/np.sqrt(npts))    

    return t_binned, y_binned, dy_binned


def prep_lc(t, y, n_std=5, detrend="wotan", wind=0.1, lim=1000, diag=False, ticid=None, cam=None,
            ccd=None, coord=None, output_dir=None, dy=None, poly_deg=10):
    if detrend != "polyfit":
        from wotan import flatten
    if diag:
        from Period_Finding import BLS

    flag = False
    
    # >> sort time array
    inds = np.argsort(t)
    t, y = t[inds], y[inds]

    # >> remove nans
    inds = np.nonzero(~np.isnan(y))
    t, y = t[inds], y[inds]
    return_dy=False
    if type(dy) != type(None):
        return_dy=True
        dy = dy[inds]

    # >> sigma-clip         
    med = np.median(y)
    std = np.std(y)
    inds = np.nonzero( (y > med - n_std*std) * (y < med + n_std*std) )
    t, y = t[inds], y[inds]
    if return_dy:
        dy = dy[inds]    

    if diag:
        if not return_dy:
            dy = np.ones(y.shape)        
        t, y, dy, period, bls_power_best, freqs, power, dur, epo = \
            BLS(t,y,dy,pmin=7,pmax=0.25,qmin=0.005,qmax=0.2,remove=True)
        prefix = 'TIC%016d'%ticid+'_0_cam_'+str(cam)+'_ccd_'+str(ccd)+\
            '_pow_'+str(bls_power_best)+'_per_'+str(round(period*1440,2))+\
            '_dur_'+str(dur)+'_epo_'+str(epo)+\
            '_ra_{}_dec_{}_'.format(coord[0], coord[1])                
        make_phase_curve(t, y, period, dy=dy, output_dir=output_dir,
                             prefix=prefix, freqs=freqs, power=power,
                             ticid=ticid, bins=100)

    if detrend == "polyfit":
        pval = np.polyfit(t, y, poly_deg)
        y_poly = np.polyval(pval, t)
        y = y - y_poly
        y += np.median(y_poly)
            
    # >> normalize    
    if np.min(y) < 0:
        y += 1 - np.min(y)
    med = np.median(y)
    y = y / med
    if return_dy:
        dy = dy/med

    # >> detrending 
    if detrend != "polyfit":
        y = flatten(t, y, window_length=wind, method='biweight')
        inds = np.nonzero(~np.isnan(y))
        t, y = t[inds], y[inds]
        if return_dy:
            dy = dy[inds]
    
    if y.shape[0] < lim:
        flag = True
    
    if return_dy:
        return t, y, dy, flag
    return t, y, flag
    

def make_phase_curve(t, y, period, dy=None, output_dir=None, prefix='', freqs=None, power=None, ticid=None,
                     bins=None, dur=None, epo=None, bls=True):
    '''Plot power spectrum and phase-folded light curve.'''
    from astropy.timeseries import TimeSeries
    from astropy.time import Time
    import astropy.units as u

    # >> normalize
    med = np.median(y)
    if np.median(y) < 0:
        y = y / np.abs(med) + 2.
        if type(dy) != type(None):
            dy = dy / np.abs(med) + 2.
    else:
        y = y/med    
        if type(dy) != type(None):
            dy = dy / med    

    # >> fold
    folded_t = t % period 
    inds = np.argsort(folded_t)
    folded_t, folded_y = folded_t[inds], y[inds]
    if type(dy) != type(None):
        folded_dy = dy[inds]
    else:
        folded_dy = None

    if type(bins) != type(None):
        folded_t, folded_y, folded_dy = bin_timeseries(folded_t, folded_y, bins, dy=folded_dy)
    
    if output_dir: # -- make plot ----------------------------------------------

        # -- initialize figure -------------------------------------------------
        fig = plt.figure(figsize=(8, 10), constrained_layout=True)
        gs = fig.add_gridspec(nrows=4, ncols=2)
        ax0_L = fig.add_subplot(gs[0, 0])
        ax0_R = fig.add_subplot(gs[0, 1])
        ax1 = fig.add_subplot(gs[1, :])
        ax2 = fig.add_subplot(gs[2, :])    
        ax3 = fig.add_subplot(gs[3, :])
        ax = [ax0_L, ax1, ax2, ax3, ax0_R]
        # fig, ax = plt.subplots(nrows=4, figsize=(8, 10))
        # ----------------------------------------------------------------------

        if type(ticid) != type(None):
            ax[0].set_title('TIC '+str(ticid)+'\nperiod: '+str(round(period*1440,2))+' min')
        else:
            ax[0].set_title('period: '+str(round(period*1440,2))+' min')
        ax[0].plot(freqs, power, '.k', ms=1, alpha=0.5)
        ax[0].set_xlabel('Frequency [1/days]')

        wind = 1000
        centr = np.argmin(np.abs(freqs - 1./period)) # >> ind of peak center
        # >> threshold power (50% of peak)
        ax[4].plot(freqs[centr-wind:centr+wind], power[centr-wind:centr+wind], '.k', ms=1)
        ax[4].set_xlabel('Frequency [1/days]') 

        thr_pow = np.median(power) + (power[centr] - np.median(power))/2
        try:
            thr_L = int(centr-wind + np.where(power[centr-wind:centr] < thr_pow)[0][-1])
            thr_R = int(centr + np.where(power[centr:centr+wind] < thr_pow)[0][0])
            ax[4].plot(freqs[thr_L:thr_R], power[thr_L:thr_R], '.r', ms=1)
        except:
            thr_R, thr_L = 0, 0

        ax[1].plot(1440/freqs, power, '.k', ms=1, alpha=0.5)
        ax[1].set_xlim([0,120])
        ax[1].set_xlabel('Period [minutes]')

        if bls:
            ax[0].set_ylabel('BLS Power')        
            ax[1].set_ylabel('BLS Power')
        else:
            ax[0].set_ylabel('LS Power')        
            ax[1].set_ylabel('LS Power')

            
        folded_t, folded_dy = np.array(folded_t), np.array(folded_dy)
        shift = np.max(folded_t) - np.min(folded_t)        
        if type(bins) != type(None):
            ax[2].errorbar(folded_t*1440, folded_y, yerr=folded_dy,
                           fmt='.k', ms=1, elinewidth=1)
            ax[2].errorbar((folded_t+shift)*1440, folded_y, yerr=folded_dy,
                           fmt='.k', ms=1, elinewidth=1)

        else:
            ax[2].plot(folded_t*1440, folded_y, '.k', ms=1)
            ax[2].plot((folded_t+shift)*1440, folded_y, '.k', ms=1)            
        ax[2].set_xlabel('Time [minutes]')
        ax[2].set_ylabel('Relative Flux')

        ax[3].plot(t, y, '.k', ms=0.8, alpha=0.8)
        ax[3].set_xlabel('Time [TJD]')
        ax[3].set_ylabel('Relative Flux')

        fig.tight_layout()
        prefix = 'wid_'+str(thr_R - thr_L)+'_'+ prefix
        plt.savefig(output_dir+prefix+'phase_curve.png', dpi=300)
        print('Saved '+output_dir+prefix+'phase_curve.png')
        plt.close()


        np.save(output_dir+prefix+'bls_power.npy',
                np.array([freqs[centr-wind:centr+wind],
                          power[centr-wind:centr+wind]]).T )
        np.save(output_dir+prefix+'phase_curve.npy',
                np.array([folded_y, folded_dy]).T)

    else:
        return np.array(folded_t), np.array(folded_y), np.array(folded_dy)

def plot_phase_curve(t, y, per, out_dir, bins=200, prefix=''):
    inds = np.nonzero(~np.isnan(y))
    t, y = t[inds], y[inds]
    folded_t, folded_y, folded_dy = make_phase_curve(t, y, per, bins=bins)
    shift = np.max(folded_t) - np.min(folded_t)
    fig, ax = plt.subplots(figsize=(5,3))
    if type(bins) != type(None):
        ax.errorbar(folded_t*1440, folded_y, yerr=folded_dy,
                       fmt='.k', ms=1, elinewidth=1)
        ax.errorbar((folded_t+shift)*1440, folded_y, yerr=folded_dy,
                       fmt='.k', ms=1, elinewidth=1)
    else:
        ax.plot(folded_t*1440, folded_y, '.k', ms=1, alpha=0.3)
        ax.plot((folded_t+shift)*1440, folded_y, '.k', ms=1, alpha=0.3) 
    ax.set_xlabel('Time [minutes]')
    ax.set_ylabel('Relative Flux')
    fig.tight_layout()
    plt.savefig(out_dir+prefix+'phase_curve.png', dpi=300)
    print('Saved '+out_dir+prefix+'phase_curve.png')

def phase(t, freq, phi0=0.):
    phi = (t * freq - phi0)
    phi -= np.floor(phi)

    return phi

def plot_orbital_phase_curve(ax, t, y, dy, freq, q, phi0, **kwargs):
    w = np.power(dy, -2)
    w /= sum(w)

    phi = phase(t, freq, phi0=phi0)
    transit = phi < q

    def ybar(mask):
        return np.dot(w[mask], y[mask]) / sum(w[mask])

    y0 = ybar(~transit)
    delta = y0 - ybar(transit)

    ax.scatter((phi[~transit] + phi0) % 1.0, y[~transit],
               c='k', s=1, alpha=0.5)
    ax.scatter((phi[transit] + phi0) % 1.0, y[transit],
               c='r', s=1, alpha=0.5)

    ax.set_xlim(0, 1)
    ax.set_xlabel('$\phi$ ($f = %.3f$)' % (freq))
    ax.set_ylabel('$y$')

def hr_digaram(gaia_tab, ra, dec, ax):

    
    from astropy.io import fits
    import astropy.units as u
    from astropy.coordinates import SkyCoord
    from astroquery.gaia import Gaia
    
    hdul = fits.open(gaia_tab)
    gid = hdul[1].data['source_id']
    gmag = hdul[1].data['phot_g_mean_mag']
    bp_rp = hdul[1].data['bp_rp']
    
    ax.plot(bp_rp, gmag, '.k', alpha=0.2, ms=0.05)
    ax.set_xlim([-0.6, 5])
    ax.set_ylim([21, 3.3])
    ax.set_xlabel('Gaia BP-RP')
    ax.set_ylabel('Absolute Magnitude (Gaia G)')

    Gaia.ROW_LIMIT = 5
    coord = SkyCoord(ra=ra, dec=dec,
                     unit=(u.degree, u.degree), frame='icrs')
    j = Gaia.cone_search_async(coord, radius=u.Quantity(3, u.arcsec))
    if len(j.get_results()['phot_g_mean_mag']) > 0 and \
       len(j.get_results()['bp_rp']) > 0:
        gmag_targ = j.get_results()['phot_g_mean_mag'][0]
        bprp_targ = j.get_results()['bp_rp'][0]
        ax.plot([bprp_targ], [gmag_targ], '^r')

def make_panel_plot(fname,sector,gaia_tab,wd_tab,tess_dir,atlas_dir,out_dir,bins=200,n_std=3,wind=0.1,pmin=400/60,pmax=0.15,qmin=0.01,qmax=0.05,ls=False):

    import pandas as pd
    import os
    from Period_Finding import BLS
    # fname example: pow_74.484985_per_162.49363_TIC0000000036085812_cam_1_ccd_1_dur_0.05_epo_0.6333333_ra_118.2292373099_dec_-6.24747586719_phase_curve.png
    # ['pow', '54.41536417029577', 'per', '113.89214', 'TIC0000000452954413', 'cam', '1', 'ccd', '4', 'ra', '122.07905627863', 'dec', '0.98336148012', 'phase', 'curve.png']

    if ls:
        ticid = int(fname.split('_')[4][3:])
        cam = fname.split('_')[6]
        ccd = fname.split('_')[8]
        per = float(fname.split('_')[3]) / 1440        
        ra = float(fname.split('_')[10])
        dec = float(fname.split('_')[12])    
        prefix = '_'.join(fname.split('_')[:13])
    else:
        # wid_14_pow_37.74368_snr_0.33903_per_8.88907_TIC0000000754984562_cam_1_ccd_1_dur_0.010638298_epo_0.8404255_ra_118.8871411275_dec_-1.78554163343_phase_curve.png
        ticid = int(fname.split('_')[8][3:])
        cam = fname.split('_')[10]
        ccd = fname.split('_')[12]
        per = float(fname.split('_')[7]) / 1440        
        ra = float(fname.split('_')[18])
        dec = float(fname.split('_')[20])    
        prefix = '_'.join(fname.split('_')[:21])
        
        # ticid = int(fname.split('_')[6][3:])
        # cam = fname.split('_')[8]
        # ccd = fname.split('_')[10]
        # per = float(fname.split('_')[5]) / 1440        
        # ra = float(fname.split('_')[16])
        # dec = float(fname.split('_')[18])    
        # prefix = '_'.join(fname.split('_')[:19])
        
        # ticid = int(fname.split('_')[4][3:])
        # cam = fname.split('_')[6]
        # ccd = fname.split('_')[8]
        # per = float(fname.split('_')[3]) / 1440        
        # ra = float(fname.split('_')[14])
        # dec = float(fname.split('_')[16])    
        # prefix = '_'.join(fname.split('_')[:17])

    fig = plt.figure(figsize=(10,6), constrained_layout=True)
    gs = fig.add_gridspec(nrows=3, ncols=2)
    ax0 = fig.add_subplot(gs[0, 1])
    ax1 = fig.add_subplot(gs[1, 1])
    ax2 = fig.add_subplot(gs[2, 1])    
    ax3 = fig.add_subplot(gs[:, 0])

    # -- tess phase curve ------------------------------------------------------

    suffix = '-{}-{}.npy'.format(cam, ccd)
    t = np.load(tess_dir+'ts'+suffix)
    ticid_list = np.load(tess_dir+'id'+suffix)
    ind = np.nonzero(ticid_list == ticid)[0][0]
    y = np.load(tess_dir+'lc'+suffix)[ind]
    t, y, flag = prep_lc(t, y, n_std=n_std, wind=wind)
    folded_t, folded_y, folded_dy = make_phase_curve(t, y, per, bins=bins)
    shift = np.max(folded_t) - np.min(folded_t)

    ax0.errorbar(folded_t*1440, folded_y, yerr=folded_dy,
                   fmt='.k', ms=1, elinewidth=1)
    ax0.errorbar((folded_t+shift)*1440, folded_y, yerr=folded_dy,
                   fmt='.k', ms=1, elinewidth=1)
    ax0.text(0.95, 0.05, str(np.round(per*1440,5))+' min', horizontalalignment='right',
             transform=ax0.transAxes)                    
    ax0.set_xlabel('Time [minutes]')
    ax0.set_ylabel('TESS Relative Flux')

    # -- atlas phase curve ------------------------------------------------------

    
    # >> load white dwarf catalog
    wd_cat  = pd.read_csv(wd_tab, header=None, sep='\s+')
    ind = np.nonzero(wd_cat[0].to_numpy() == ticid)[0][0]
    if not np.isnan(wd_cat.iloc[ind][3]):
        gid = int(wd_cat.iloc[ind][3])
    else:
        gid = None
    if os.path.exists(atlas_dir+str(gid)):
        data = atlas_lc(atlas_dir+str(gid))
        t, y, dy = data[:,0], data[:,1], data[:,2]
        try:
            _, _, _, period, bls_power_best, freqs, power, dur, epo, delta = \
                BLS(t,y,dy,pmin=pmin,pmax=pmax,qmin=qmin,qmax=qmax,remove=True)
        except:
            dy=np.ones(y.shape)
            _, _, _, period, bls_power_best, freqs, power, dur, epo, delta = \
                BLS(t,y,dy,pmin=pmin,pmax=pmax,qmin=qmin,qmax=qmax,remove=True)            
        prefix1 = 'ATLAS_'+prefix+'_'
        make_phase_curve(t, y, period, dy=dy, output_dir=out_dir,
                         prefix=prefix1, freqs=freqs, power=power,
                         ticid=ticid, bins=100)
        
        folded_t, folded_y, folded_dy = make_phase_curve(t, y, period, bins=bins, dy=dy)
        shift = np.max(folded_t) - np.min(folded_t)
        ax1.errorbar(folded_t*1440, folded_y, yerr=folded_dy,
                       fmt='.k', ms=1, elinewidth=1)
        ax1.errorbar((folded_t+shift)*1440, folded_y, yerr=folded_dy,
                       fmt='.k', ms=1, elinewidth=1)
        ax1.text(0.95, 0.05, str(np.round(period*1440,5))+' min', horizontalalignment='right',
                 transform=ax1.transAxes)                
    ax1.set_xlabel('Time [minutes]')
    ax1.set_ylabel('ATLAS Relative Flux')

    # -- ztf phase curve ------------------------------------------------------

    t, y, dy = [], [], []
    os.system('python /data/ZTF_Lightcurves/get_LC.py {} {} g'.format(ra, dec))
    if os.path.exists('/home/echickle/{}_{}_g.lc'.format(ra,dec)) and \
       os.path.getsize('/home/echickle/{}_{}_g.lc'.format(ra,dec)) != 0:
        data = np.loadtxt('/home/echickle/{}_{}_g.lc'.format(ra,dec))
        med = np.median(data[:,1])
        t.extend(data[:,0])
        y.extend(data[:,1] - np.median(data[:,1]))
        dy.extend(data[:,2])
    os.system('python /data/ZTF_Lightcurves/get_LC.py {} {} r'.format(ra, dec))
    if os.path.exists('/home/echickle/{}_{}_r.lc'.format(ra,dec)) and \
       os.path.getsize('/home/echickle/{}_{}_r.lc'.format(ra,dec)) != 0:    
        data  = np.loadtxt('/home/echickle/{}_{}_r.lc'.format(ra,dec))
        med = np.median(data[:,1])        
        t.extend(data[:,0])
        y.extend(data[:,1] - np.median(data[:,1]))
        dy.extend(data[:,2])
    os.system('python /data/ZTF_Lightcurves/get_LC.py {} {} i'.format(ra, dec))
    if os.path.exists('/home/echickle/{}_{}_i.lc'.format(ra,dec)) and \
       os.path.getsize('/home/echickle/{}_{}_i.lc'.format(ra,dec)) != 0: 
        data = np.loadtxt('/home/echickle/{}_{}_i.lc'.format(ra,dec))
        med = np.median(data[:,1])        
        t.extend(data[:,0])
        y.extend(data[:,1] - np.median(data[:,1]))
        dy.extend(data[:,2])
        
    if len(t) != 0:
        t, y, dy = np.array(t), np.array(y) + med, np.array(dy)
        _, _, _, period, bls_power_best, freqs, power, dur, epo, delta = \
            BLS(t,y,dy,pmin=pmin,pmax=pmax,qmin=qmin,qmax=qmax,remove=True)
        prefix1 = 'ZTF_'+prefix+'_'
        make_phase_curve(t, y, period, dy=dy, output_dir=out_dir,
                         prefix=prefix1, freqs=freqs, power=power,
                         ticid=ticid)
        folded_t, folded_y, folded_dy = make_phase_curve(t, y, period, dy=dy)
        shift = np.max(folded_t) - np.min(folded_t)
        ax2.errorbar(folded_t*1440, folded_y, yerr=folded_dy,
                       fmt='.k', ms=1, elinewidth=1)
        ax2.errorbar((folded_t+shift)*1440, folded_y, yerr=folded_dy,
                       fmt='.k', ms=1, elinewidth=1)

        ax2.text(0.95, 0.05, str(np.round(period*1440,5))+' min', horizontalalignment='right',
                 transform=ax2.transAxes)        
    ax2.set_xlabel('Time [minutes]')
    ax2.set_ylabel('ZTF Relative Flux')

            

    # -- hr diagram ------------------------------------------------------------

    hr_digaram(gaia_tab, ra, dec, ax3)
    
    fig.tight_layout()
    plt.savefig(out_dir+prefix+'_panel.png', dpi=300)
    print('Saved '+out_dir+prefix+'_panel.png')
    plt.close()


    
