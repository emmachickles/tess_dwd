import numpy as np
import matplotlib.pyplot as plt
import pdb

def load_atlas_lc(f, n_std=2, clip=True):

    ###MJD     m   dm uJy duJy F err chi/N  RA   Dec    x    y  maj min phi apfit Sky  ZP  Obs  GaiaID
    data=np.loadtxt(f,usecols=(0,3,4,16),skiprows=0)    
    Filter=np.loadtxt(f,usecols=(5),skiprows=0,dtype=str) 

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

    # >> sigma clip
    # std = np.std(y)
    # inds = np.nonzero( (y > med - n_std*std) * (y < med + n_std*std) )
    # t, y, dy = t[inds], y[inds], dy[inds] 

    if clip:
        q3, q1 = np.percentile(y, [75 ,25])
        iqr=(q3-q1)/2

        good_idx=(y-np.median(y))<3*iqr # !! 3 
        t=t[good_idx]
        dy=dy[good_idx]
        y=y[good_idx]

        good_idx=(np.median(y)-y)<10*iqr
        t=t[good_idx]
        dy=dy[good_idx]
        y=y[good_idx]

        good_idx=dy>0
        t=t[good_idx]
        y=y[good_idx]
        dy=dy[good_idx]  

    # >> convert to BJD
    coords=np.loadtxt(f,usecols=(8,9),skiprows=0)
    ra, dec = np.mean(coords[:,0]), np.mean(coords[:,1])
    t = BJDConvert(t,ra,dec, date_format='mjd').value    

    return t, y, dy, ra, dec

def get_atlas_lc(ticid, wd_tab, atlas_dir):
    import pandas as pd
    import os
    # >> load white dwarf catalog
    wd_cat  = pd.read_csv(wd_tab, header=None, sep='\s+', dtype='str')

    # >> find gaia id 
    ind = np.nonzero(wd_cat[0].to_numpy().astype('int') == ticid)[0][0]
    if type(wd_cat.iloc[ind][3]) == type(""):
        gid = str(wd_cat.iloc[ind][3])
    else:
        gid = None
    if os.path.exists(atlas_dir+str(gid)):
        return atlas_dir+str(gid)
    else:
        return None

def get_ztf_lc(ra, dec):
    import os
    os.system('python /data/ZTF_Lightcurves/get_LC.py {} {} g'.format(ra, dec))    
    os.system('python /data/ZTF_Lightcurves/get_LC.py {} {} r'.format(ra, dec))
    os.system('python /data/ZTF_Lightcurves/get_LC.py {} {} i'.format(ra, dec))

    fnames = []
    fname_g = '/home/echickle/{}_{}_g.lc'.format(ra,dec)
    if os.path.exists(fname_g):
        if os.path.getsize(fname_g) != 0:
            fnames.append(fname_g)
    fname_r = '/home/echickle/{}_{}_r.lc'.format(ra,dec)
    if os.path.exists(fname_r):
        if os.path.getsize(fname_r) != 0:
            fnames.append(fname_r)
    fname_i = '/home/echickle/{}_{}_i.lc'.format(ra,dec)
    if os.path.exists(fname_i):
        if os.path.getsize(fname_i) != 0:
            fnames.append(fname_i)

    return fnames

def load_ztf_lc(fnames, n_std=5):
    t, y, dy = [], [], []
    
    for i in range(len(fnames)):
        f = fnames[i]
        data = np.loadtxt(f)
        med = np.median(data[:,1])
        t.extend(data[:,0])
        y.extend(data[:,1] - np.median(data[:,1]))
        dy.extend(data[:,2])    

    t, y, dy = np.array(t), np.array(y) + med, np.array(dy)

    # std = np.std(y)
    # inds = np.nonzero( (y > med - n_std*std) * (y < med + n_std*std) )
    # t, y, dy = t[inds], y[inds], dy[inds]  

    q3, q1 = np.percentile(y, [75 ,25])
    iqr=(q3-q1)/2

    good_idx=(y-np.median(y))<3*iqr
    t=t[good_idx]
    dy=dy[good_idx]
    y=y[good_idx]

    good_idx=(np.median(y)-y)<10*iqr
    t=t[good_idx]
    dy=dy[good_idx]
    y=y[good_idx]

    return t, y, dy
        
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
    if len(t) < bins*2:
        t_binned, y_binned = t, y
        if dy is None:
            dy_binned = np.zeros(y.shape)
        else:
            dy_binned = dy
    else:    
        # >> should optimize this with scipy.stats.binned_statistic
        inds = np.argsort(t)
        t, y = t[inds], y[inds]
        if dy is not None:
            dy = dy[inds]

        trunc = len(y) // bins * bins
        t_binned = np.split(t[:trunc], bins)
        y_binned = np.split(y[:trunc], bins)
        if dy is not None:
            dy_binned = np.split(dy[:trunc], bins)

        if trunc < len(y):
            y_binned[-1] = np.append(y_binned[-1], y[trunc:])
            t_binned[-1] = np.append(t_binned[-1], t[trunc:])
            if dy is not None:
                dy_binned[-1] = np.append(dy_binned[-1], dy[trunc:])

        if dy is None:
            dy_binned = []
        for i in range(bins):
            t_binned[i] = np.average(t_binned[i])

            npts = len(y_binned[i])

            if dy is None:
                err = np.std(y_binned[i]) / np.sqrt(npts)            
                y_binned[i] = np.average(y_binned[i])
                dy_binned.append(err)
            else:
                y_binned[i] = np.average(y_binned[i], weights=1./dy_binned[i]**2)
                dy_binned[i] = np.average(dy_binned[i]/np.sqrt(npts))

    return np.array(t_binned), np.array(y_binned), np.array(dy_binned)

def normalize_lc(y, dy=None):
    if np.min(y) < 0:
        y -= np.min(y)
    med = np.median(y)
    y = y/med
    if dy is not None:
        dy = dy/med
    return y, dy

def prep_lc(t, y, n_std=5, detrend="wotan", wind=0.1, lim=1000, ticid=None, cam=None,
            ccd=None, coord=None, output_dir=None, dy=None, poly_deg=10):
    if detrend != "polyfit":
        from wotan import flatten

    flag = False
    
    # >> sort time array
    inds = np.argsort(t)
    t, y = t[inds], y[inds]

    # >> remove nans
    inds = np.nonzero(~np.isnan(y))
    t, y = t[inds], y[inds]

    if dy is not None:
        dy = dy[inds]

    # >> sigma-clip         
    med = np.median(y)
    std = np.std(y)
    inds = np.nonzero( (y > med - n_std*std) * (y < med + n_std*std) )
    t, y = t[inds], y[inds]
    if dy is not None:
        dy = dy[inds]    

    if detrend == "polyfit":
        pval = np.polyfit(t, y, poly_deg)
        y_poly = np.polyval(pval, t)
        y = y - y_poly
        y += np.median(y_poly)

    if y.shape[0] < lim:
        flag = True
        if dy is not None:
            return t, y, dy, flag
        return t, y, flag
            
    # >> normalize
    y, dy = normalize_lc(y, dy)

    # >> detrending 
    if detrend != "polyfit":
        y = flatten(t, y, window_length=wind, method='biweight')
        inds = np.nonzero(~np.isnan(y))
        t, y = t[inds], y[inds]
        if dy is not None:
            dy = dy[inds]
    
    if y.shape[0] < lim:
        flag = True

    if dy is not None:
        return t, y, dy, flag
    return t, y, flag
    

def calc_snr(t, y, period, q, phi0):
    from scipy.stats import median_abs_deviation

    t_mean = np.mean(t)
    t = t - t_mean

    # >> get epoch (center of first eclipse in time series)
    epo = (phi0 + 0.5 * q) * period
    epo += period + int((np.min(t) - epo) / period)*period
    epo_TJD = epo + t_mean

    # -- phase fold ------------------------------------------------------------
    phi= np.mod(t - epo, period) / period # >> phase
    phi[phi > 0.5] -= 1 # >> center transit at phase=0    

    # -- calculate SNR ---------------------------------------------------------
    # >> in-eclipse datapoints:
    transit = np.abs(phi) < 0.5*q
    # >> in-eclipse and nearby datapoints:
    near_transit = np.abs(phi) < 1.5*q
    # >> other out-of-eclipse datapoints:
    out_transit = np.abs(phi) > q
 
    # >> get delta (depth of transit)
    # avg_out = np.mean(y[out_transit])
    # avg_in = np.mean(y[transit])
    # err_out = np.std(y[out_transit]) / np.sqrt(len(y[out_transit]))
    # err_in = np.std(y[transit]) / np.sqrt(len(y[transit]))
    # err = np.sqrt(err_out**2 + err_in**2)
    
    if np.count_nonzero(transit) > 0:
        avg_out = np.median(y[out_transit])
        avg_in = np.median(y[transit])
        # err_out = median_abs_deviation(y[out_transit])
        # err_in = median_abs_deviation(y[transit])
        # err = np.sqrt(err_out**2 + err_in**2)
        err = median_abs_deviation(y[out_transit])

        delta = np.abs(avg_out - avg_in)
        snr = delta/err
    else:
        snr = 0.

    nt = np.count_nonzero(transit)
    dphi = np.max( np.diff( np.sort(phi) ) )
    return snr, phi, transit, near_transit, epo_TJD, nt, dphi

def calc_peak_stats(freqs, power, nearpeak=3000):
    from scipy.stats import median_abs_deviation
    peak = np.argmax(power)
    
    # -- significance of peak ----------------------------------------------
    sig=(np.max(power)-np.median(power)) / median_abs_deviation(power)

    # -- width of peak -----------------------------------------------------
    # >> threshold power (50% of peak)
    med_near = np.median(power[max(0,peak-nearpeak) : peak+nearpeak])
    thr_pow = med_near + (power[peak] - med_near)/2
    inds_L = np.where(power[max(0,peak-nearpeak):peak] < thr_pow)[0]
    if len(inds_L) > 0:
        thr_L = int(max(0,peak-nearpeak) + inds_L[-1])
    else:
        thr_L = peak
    inds_R = np.where(power[peak:peak+nearpeak] < thr_pow)[0]
    if len(inds_R) > 0:
        thr_R = int(peak + inds_R[0])
    else:
        thr_R = peak

    wid = thr_R - thr_L

    f_best = freqs[peak]
    period=1.0/f_best

    return peak, sig, wid, f_best, period

def calc_sine_fit(t, y, period):
    from scipy.optimize import curve_fit
    t_folded = t%period
    w = 2*np.pi / period
    def sinfunc(t, A, p, c):  return A * np.sin(w*t + p) + c
    popt, pcov = curve_fit(sinfunc, t_folded, y)
    A, p, c = popt # >> amplitude, phase, offset
    fitfunc = sinfunc(t_folded, A, p, c)
    err = np.mean( ( fitfunc - y ) ** 2)
    return err, fitfunc
    
def vet_plot(t, y, freqs, power, q=None, phi0=None, dy=None, output_dir=None, suffix='',
             objid=None, objid_type='TICID', bins=100, bls=True, save_npy=False, nearpeak=3000,
             wd_tab='WDs.txt', wd_main='GaiaEDR3_WD_main.fits', rp_ext='GaiaEDR3_WD_RPM_ext.fits',
             snr_threshold=0, pow_threshold=0, per_threshold=14400, wid_threshold=0):
    '''Plot power spectrum and phase-folded light curve.
    * q : ratio of eclipse duration to period
    * phi0 : start of eclipse phase
    * nearpeak: frequency bins near peak'''
    from astropy.timeseries import TimeSeries
    from astropy.time import Time
    import astropy.units as u

    # == normalize =============================================================
    y, dy = normalize_lc(y, dy)

    # == vetting ===============================================================

    # -- peak statistics -------------------------------------------------------
    peak, sig, wid, f_best, period = calc_peak_stats(freqs, power, nearpeak=nearpeak)

    # -- calculate SNR ---------------------------------------------------------
    if bls:
        snr, phi, transit, near_transit, epo, nt, dphi = calc_snr(t, y, period, q, phi0)
        
    # -- calculate fit ---------------------------------------------------------
    if not bls:
        fit, fitfunc = calc_sine_fit(t, y, period)
        
    # -- calculate companion radius --------------------------------------------
    if bls:
        # >> semi-major axis in km, from Kepler's third law
        a = (0.6 * (period / 365.25)**2)**(1/3) * 1.496e8

        dur = q*period*1440 # >> duration in minutes
        vel = 2*np.pi*a / (period * 86400) # >> velocity in km/s
        rp = vel*dur*60/2 # >> radius in km
        rp = rp / 6370 # >> radius in Earth radii

    # == bin ===================================================================
    folded_t, folded_y, folded_dy = bin_timeseries(t%period, y, bins, dy=dy)
    
    # print('Fold light curve done!')
    # == make plot =============================================================
    if sig > pow_threshold and snr > snr_threshold and period*1440 < per_threshold and\
       wid > wid_threshold: 

        # -- initialize figure -------------------------------------------------
        plot_pg=True
        if len(freqs) < 1e6:
            plot_pg=False

        if plot_pg:
            fig = plt.figure(figsize=(8, 10), constrained_layout=True)
            gs = fig.add_gridspec(nrows=4, ncols=2)
            ax0_L = fig.add_subplot(gs[0, 0])
            ax0_R = fig.add_subplot(gs[0, 1])
            ax1 = fig.add_subplot(gs[1, :])            
            ax2 = fig.add_subplot(gs[2, :])    
            if bls:
                ax3_L = fig.add_subplot(gs[3, 0])
                ax3_R = fig.add_subplot(gs[3, 1])
            else:
                ax3 = fig.add_subplot(gs[3, :])            
        else:
            fig = plt.figure(figsize=(8, 8), constrained_layout=True)
            gs = fig.add_gridspec(nrows=3, ncols=2)
            ax0_L = fig.add_subplot(gs[0, 0])
            ax0_R = fig.add_subplot(gs[0, 1])
            ax2 = fig.add_subplot(gs[1, :])    
            if bls:
                ax3_L = fig.add_subplot(gs[2, 0])
                ax3_R = fig.add_subplot(gs[2, 1])
            else:
                ax3 = fig.add_subplot(gs[2, :])

        # -- title -------------------------------------------------------------
        if bls:
            suptext ='period: '+str(round(period*1440,2))+' min\nduration: '+\
                      str(round(dur,2))+' mins, epo: '+\
                      str(round(epo,5))+' TJD, snr: '+\
                      str(round(snr, 2))+'\nradius assuming .6 '+r'$M_\odot$'+\
                      ' WD: '+str(round(rp, 2))+' '+r'$R_\oplus$'
        else:
            suptext='period: '+str(round(period*1440,2))+' min'            
        if objid is not None:
            suptext = objid_type+' '+str(objid)+', ' + suptext
        fig.suptitle(suptext)

        # -- plot periodogram --------------------------------------------------
        # if len(freqs) < 1e6:
        #     ax0_L.plot(freqs, power, '.k', ms=1, alpha=0.5, rasterized=True)
        #     ax0_L.set_xlabel('Frequency [1/days]')
        hr_diagram_wd(objid, objid_type, ax0_L, wd_tab=wd_tab, wd_main=wd_main,
                      rp_ext=rp_ext)

        # >> threshold power (50% of peak)
        ax0_R.plot(freqs[max(0,peak-nearpeak):peak+nearpeak],
                   power[max(0,peak-nearpeak):peak+nearpeak], '.k', ms=1)
        ax0_R.plot(freqs[peak-wid//2:peak+wid//2],
                   power[peak-wid//2:peak+wid//2], '.r', ms=1)
        ax0_R.set_xlabel('Frequency [1/days]')

        if plot_pg:
            ax1.plot(1440/freqs, power, '.k', ms=1, alpha=0.5, rasterized=True)
            ax1.set_xlim([np.min(1440/freqs), np.max(1440/freqs)])
            ax1.set_xlabel('Period [minutes]')

            if bls:
                # ax0_L.set_ylabel('BLS Power')        
                ax1.set_ylabel('BLS Power')
            else:
                # ax0_L.set_ylabel('LS Power')        
                ax1.set_ylabel('LS Power')

        # -- plot phase curve --------------------------------------------------


        folded_t, folded_dy = np.array(folded_t), np.array(folded_dy)
        shift = np.max(folded_t) - np.min(folded_t)        
        if type(bins) != type(None):
            ax2.errorbar(folded_t*1440, folded_y, yerr=folded_dy,
                           fmt='.k', ms=1, elinewidth=1)
            ax2.errorbar((folded_t+shift)*1440, folded_y, yerr=folded_dy,
                           fmt='.k', ms=1, elinewidth=1)

        else:
            ax2.plot(folded_t*1440, folded_y, '.k', ms=1)
            ax2.plot((folded_t+shift)*1440, folded_y, '.k', ms=1)            
        ax2.set_xlabel('Time [minutes]')
        ax2.set_ylabel('Relative Flux')

        # -- plot full light curve ---------------------------------------------
        if bls:
            ax3_L.plot(t, y, '.k', ms=0.8, alpha=0.8)   
            ax3_L.set_xlabel('Time [TJD]')
            ax3_L.set_ylabel('Relative Flux')

            
            # >> plot phase curve zoomed in on transit 
            ax3_R.axvline(-0.5*q*period*1440, color='k', lw=0.5, ls='dashed')
            ax3_R.axvline(0.5*q*period*1440, color='k', lw=0.5, ls='dashed')
            t_transit = phi[near_transit]*period*1440

            ax3_R.plot(t_transit, y[near_transit], '.k', ms=1)
            ax3_R.set_xlim([np.min(t_transit), np.max(t_transit)])

            w = max(1, int(0.1*np.count_nonzero(near_transit)))
            inds = np.argsort(phi[near_transit])
            if np.count_nonzero(near_transit) > 0:
                phiconv = np.convolve(phi[near_transit][inds], np.ones(w), 'valid') / w
                yconv = np.convolve(y[near_transit][inds], np.ones(w), 'valid') / w
                ax3_R.plot(phiconv*period*1440, yconv, '-')
                # ax3_R.set_ylim([np.min(yconv)-0.1, np.max(yconv)+0.1])

            ax3_R.set_xlabel('Time [minutes]')
            ax3_R.set_ylabel('Relative Flux')
        else:
            ax3.plot(t, y, '.k', ms=0.8, alpha=0.8)
            ax3.set_xlabel('Time [TJD]')
            ax3.set_ylabel('Relative Flux')

        # -- save figure -------------------------------------------------------
        fig.tight_layout()

        if bls:
            prefix = 'pow_'+str(round(sig, 2))+'_snr_'+str(round(snr,2))+'_wid_'+\
                     str(wid)+'_per_'+str(round(period*1440,8))+'_q_'+\
                     str(round(q,5))+\
                     '_phi0_'+str(round(phi0,5))
        else:
            prefix = 'pow_'+str(round(sig, 2))+'_wid_'+str(wid)+'_per_'+\
                     str(round(period*1440,8))

        if bls:
            plt.savefig(output_dir+prefix+suffix+'_bls.png', dpi=100)
            print('Saved '+output_dir+prefix+suffix+'_bls.png')
        else:
            plt.savefig(output_dir+prefix+suffix+'_ls.png', dpi=300)
            print('Saved '+output_dir+prefix+suffix+'_ls.png')
        plt.close()


        if save_npy:
            np.save(output_dir+prefix+'_bls_power.npy',
                    np.array([freqs[max(0,peak-nearpeak):peak+nearpeak],
                              power[max(0,peak-nearpeak):peak+nearpeak]]).T )
            np.save(output_dir+prefix+'_phase_curve.npy',
                    np.array([folded_y, folded_dy]).T)

    # -- return values ---------------------------------------------------------
    if bls:
        return sig, snr, wid, period, period*1440, q, phi0, dur, epo, rp, nt, dphi
    else:
        return sig, wid, period, period*1440

            
def plot_phase_curve(ax, folded_t, folded_y, folded_dy, period=None,
                     ylabel="Relative Flux"):
    shift = np.max(folded_t) - np.min(folded_t)    
    ax.errorbar(folded_t*1440, folded_y, yerr=folded_dy, fmt=".k", ms=1,
                elinewidth=1)
    ax.errorbar((folded_t+shift)*1440, folded_y, yerr=folded_dy, fmt=".k", ms=1,
                elinewidth=1)
    if period is not None:
        ax.text(0.95, 0.05, str(np.round(period*1440,5))+" min",
                horizontalalignment="right", transform=ax.transAxes)                    
    ax.set_xlabel("Time [minutes]")
    if ylabel is not None:
        ax.set_ylabel(ylabel)

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

def BJDConvert(times, RA, Dec, date_format='mjd', telescope='Palomar'):
    '''Function for converting a series of timestamps to Barycentric Julian
    Date format in Barycentric Dynamical time'''
    import numpy as np
    from astropy.time import Time
    from astropy.coordinates import EarthLocation
    from astropy.coordinates import SkyCoord  # High-level coordinates
    from astropy.coordinates import ICRS, Galactic, FK4, FK5, BarycentricTrueEcliptic  # Low-level frames
    from astropy.coordinates import Angle, Latitude, Longitude  # Angles
    import astropy.units as u
    
    t = Time(times,format=date_format,scale='utc')
    t2=t.tcb
    c = SkyCoord(RA,Dec, unit="deg")
    d=c.transform_to(BarycentricTrueEcliptic)
    Observatory=EarthLocation.of_site(telescope)
    delta=t2.light_travel_time(c,kind='barycentric',location=Observatory)
    BJD_TCB=t2+delta
    return BJD_TCB

def hr_diagram_wd(objid, objid_type, ax, wd_tab='WDs.txt', wd_main='/data/GaiaEDR3_WD_main.fits',
                  rp_ext='/data/GaiaEDR3_WD_RPM_ext.fits'):
    '''Plot white dwarf track and surrounding areas.'''
    from astropy.io import fits
    import time

    if objid_type=='GAIAID':
        gid = np.int64(objid)
    elif objid_type=='TICID':
        wd_cat  = pd.read_csv(wd_tab, header=None, sep='\s+', dtype='str')
        ind = np.nonzero(wd_cat[0].to_numpy().astype('int') == objid)[0][0]
        # if type(wd_cat.iloc[ind][3]) == type(""):
        gid = np.int64(wd_cat.iloc[ind][3])
    
    source_id = np.empty(0)
    bp_rp = np.empty(0)
    parallax = np.empty(0)
    gmag = np.empty(0)

    # start =time.time()
    maincat = fits.open(wd_main)
    source_id = np.append(source_id, maincat[1].data['source_id'])
    bp_rp = np.append(bp_rp, maincat[1].data['bp_rp'])
    parallax = np.append(parallax, maincat[1].data['parallax'])
    gmag = np.append(gmag, maincat[1].data['phot_g_mean_mag'])
    # end=time.time()
    # print(end-start)

    # start = time.time()
    rpmext = fits.open(rp_ext) 
    source_id = np.append(source_id, rpmext[1].data['source_id'])
    bp_rp = np.append(bp_rp, rpmext[1].data['bp_rp'])
    parallax = np.append(parallax, rpmext[1].data['parallax'])
    gmag = np.append(gmag, rpmext[1].data['phot_g_mean_mag'])
    # end = time.time()
    # print(end-start)

    # start=time.time()
    abs_mag = gmag+5*(np.log10(parallax)-2)
    # end = time.time()
    # print(end-start)

    ind = np.nonzero(source_id == gid)[0][0]

    # ax.plot(bp_rp, abs_mag, '.k', alpha=0.2, ms=0.05)
    # ax.set_xlim([-0.6, 1.9])
    # ax.set_ylim([15.5, 4])

    # start=time.time()
    _ = ax.hist2d(bp_rp, abs_mag, bins=1000, range=[[-0.6, 1.9], [4, 15.5]], density=True,
                  cmin=0.03)
    ax.set_xlabel('Gaia BP-RP')
    ax.set_ylabel('Absolute Magnitude (Gaia G)')
    ax.invert_yaxis()
    if not np.isnan(abs_mag[ind]):
        ax.plot([bp_rp[ind]], [abs_mag[ind]], '^r')
        ax.text(0.95, 0.95, "bp_rp: "+str(round(bp_rp[ind],2))+\
                "\ng_mean_mag: "+str(round(gmag[ind], 2))+\
                "\nparallax: "+str(round(parallax[ind],2))+\
                "\nabs_mag: "+str(round(abs_mag[ind], 2)),
                ha="right", va='top',transform=ax.transAxes)                 
    else:
        ax.text(0.95, 0.95, "bp_rp: "+str(bp_rp[ind])+\
                "\ng_mean_mag: "+str(gmag[ind])+\
                "\nparallax: "+str(parallax[ind]),
                ha="right", va='top', transform=ax.transAxes)

    # end=time.time()
    # print(end-start)
    
    
def hr_diagram(gaia_tab, ra, dec, ax):

    
    from astropy.io import fits
    import astropy.units as u
    from astropy.coordinates import SkyCoord
    from astroquery.gaia import Gaia
    
    hdul = fits.open(gaia_tab)
    gid = hdul[1].data['source_id']
    bp_rp = hdul[1].data['bp_rp']

    parallax = hdul[1].data['parallax']
    gmag = hdul[1].data['phot_g_mean_mag']
    abs_mag = gmag+5*(np.log10(parallax)-2)
    
    ax.plot(bp_rp, abs_mag, '.k', alpha=0.2, ms=0.05)
    ax.set_xlim([-0.6, 5])
    ax.set_ylim([17.5, 0])
    ax.set_xlabel('Gaia BP-RP')
    ax.set_ylabel('Absolute Magnitude (Gaia G)')

    Gaia.ROW_LIMIT = 5
    Gaia.MAIN_GAIA_TABLE = "gaiadr3.gaia_source"
    coord = SkyCoord(ra=ra, dec=dec,
                     unit=(u.degree, u.degree), frame='icrs')
    j = Gaia.cone_search_async(coord, radius=u.Quantity(3, u.arcsec))
    if len(j.get_results()['phot_g_mean_mag']) > 0 and \
       len(j.get_results()['bp_rp']) > 0:
        bprp_targ = j.get_results()['bp_rp'][0]        
        apparent_mag = j.get_results()['phot_g_mean_mag'][0]
        parallax = j.get_results()['parallax'][0]
        if str(parallax) == '--':
            abs_mag = None
            ax.text(0.95, 0.05, "bp_rp: "+str(bprp_targ)+\
                    "\ng_mean_mag: "+str(apparent_mag)+\
                    "\nparallax: "+str(parallax),
                    horizontalalignment="right", transform=ax.transAxes)
        else:
            abs_mag = apparent_mag+5*(np.log10(parallax)-2)        
            ax.plot([bprp_targ], [abs_mag], '^r')

            ax.text(0.95, 0.05, "bp_rp: "+str(round(bprp_targ,2))+\
                    "\ng_mean_mag: "+str(round(apparent_mag, 2))+\
                    "\nparallax: "+str(parallax)+\
                    "\nabs_mag: "+str(round(abs_mag, 2)),
                    horizontalalignment="right", transform=ax.transAxes) 

def make_panel_plot(fname_atlas,fnames_ztf,tess_dir,ticid,cam,ccd,
                    per,ra,dec,gaia_tab,out_dir,suffix,
                    per_tess=None,per_atlas=None,per_ztf=None,
                    bins=100,n_std=5,wind=0.1,qmin=0.01,
                    qmax=0.15,bls=False):

    import pandas as pd
    import os
    from Period_Finding import BLS
        
    fig = plt.figure(figsize=(16,6), constrained_layout=True)
    gs = fig.add_gridspec(nrows=3, ncols=4)
    ax0 = fig.add_subplot(gs[0, 1])
    ax1 = fig.add_subplot(gs[1, 1])
    ax2 = fig.add_subplot(gs[2, 1])    
    ax3 = fig.add_subplot(gs[:, 0])
    if per_tess is None:
        ax4 = fig.add_subplot(gs[0, 2])
        ax7 = fig.add_subplot(gs[0, 3])
    if per_atlas is None:
        ax5 = fig.add_subplot(gs[1, 2])
        ax8 = fig.add_subplot(gs[1, 3])
    if per_ztf is None:
        ax6 = fig.add_subplot(gs[2, 2])
        ax9 = fig.add_subplot(gs[2, 3])

    # -- tess phase curve ------------------------------------------------------

    f_suffix = '-{}-{}.npy'.format(cam, ccd)
    t = np.load(tess_dir+'ts'+f_suffix)
    ticid_list = np.load(tess_dir+'id'+f_suffix)
    ind = np.nonzero(ticid_list == ticid)[0][0]
    y = np.load(tess_dir+'lc'+f_suffix)[ind]
    t, y, flag = prep_lc(t, y, n_std=n_std, wind=wind)
    
    if per_tess is None: # -- tess period finding ------------------------------
        dy = np.ones(y.shape) * 0.1
        if bls:
            t, y, dy, per_tess, bls_power_best, freqs, power, q, phi0 = \
                BLS(t,y,dy,pmin=400/60,pmax=10,qmin=qmin,qmax=qmax,remove=True)
        else:
            _, _, _, per_tess, ls_power_best, freqs, power = \
                LS_Astropy(t,y,dy,pmax=10)        
        prefix1 = 'TESS_'
        vet_plot(t, y, freqs, power, q, phi0, output_dir=out_dir+prefix1,
                 objid=ticid, suffix='_'+suffix)

        ax4.plot(1440/freqs, power, '.k', ms=1, alpha=0.5)
        ax4.set_xlabel('Period [minutes]')

        centr = np.argmin(np.abs(freqs - 1/per))
        ax7.axvline(x=1/per, color='r', linestyle='dashed')
        ax7.plot(freqs[max(0,centr-3000):centr+3000],
                 power[max(0,centr-3000):centr+3000], '.k', ms=1)
        ax7.set_xlabel('Frequency [1/days]')

        if bls:
            ax4.set_ylabel('TESS BLS Power')        
            ax7.set_ylabel('TESS BLS Power')
        else:
            ax4.set_ylabel('TESS LS Power') 
            ax7.set_ylabel('TESS LS Power')

        
    y, _ = normalize_lc(y)
    folded_t, folded_y, folded_dy = bin_timeseries(t%per_tess, y, bins)
    plot_phase_curve(ax0, folded_t, folded_y, folded_dy, period=per_tess,
                     ylabel="TESS Relative Flux")

        
    # -- atlas phase curve ------------------------------------------------------

    if type(fname_atlas) == type(""):
        t, y, dy, _, _ = load_atlas_lc(fname_atlas)
        if per_atlas is None: # -- atlas period finding ------------------------
            if bls:
                t, y, dy, per_atlas, bls_power_best, freqs, power, q, phi0 = \
                    BLS(t,y,dy,pmin=2,pmax=10,qmin=qmin,qmax=qmax,remove=False)
            else:
                _, _, _, per_atlas, ls_power_best, freqs, power = \
                    LS_Astropy(t,y,dy,pmax=10)     

            prefix1 = 'ATLAS_'
            vet_plot(t, y, freqs, power, q, phi0, output_dir=out_dir+prefix1,
                     objid=ticid, dy=dy, suffix='_'+suffix)

            ax5.plot(1440/freqs, power, '.k', ms=1, alpha=0.5)
            ax5.set_xlabel('Period [minutes]')

            ax8.axvline(x=1/per, color='r', linestyle='dashed')        
            ax8.plot(freqs, power, '.k', ms=1)
            ax8.set_xlabel('Frequency [1/days]')
            ax8.set_xlim(ax7.get_xlim())

            if bls:
                ax5.set_ylabel('ATLAS BLS Power')        
                ax8.set_ylabel('ATLAS BLS Power')
            else:
                ax5.set_ylabel('ATLAS LS Power') 
                ax8.set_ylabel('ATLAS LS Power')            
            
        folded_t, folded_y, folded_dy = bin_timeseries(t%per_atlas, y, bins, dy=dy)
        plot_phase_curve(ax1, folded_t, folded_y, folded_dy, period=per_atlas,
                         ylabel="ATLAS Relative Flux")

    # -- ztf phase curve -------------------------------------------------------

    if len(fnames_ztf) > 0:
        t, y, dy = load_ztf_lc(fnames_ztf)
        if per_ztf is None: # -- ztf period finding ----------------------------
            if bls:
                t, y, dy, per_ztf, bls_power_best, freqs, power, q, phi0 = \
                    BLS(t,y,dy,pmin=2,pmax=10,qmin=qmin,qmax=qmax,remove=False)
            else:
                _, _, _, per_ztf, ls_power_best, freqs, power = \
                    LS_Astropy(t,y,dy,pmax=pmax)                    
            prefix1 = 'ZTF_'
            vet_plot(t, y, freqs, power, q, phi0, output_dir=out_dir+prefix1,
                     objid=ticid, dy=dy, plot_threshold=0, suffix='_'+suffix)

            ax6.plot(1440/freqs, power, '.k', ms=1, alpha=0.5)
            ax6.set_xlabel('Period [minutes]')

            ax9.axvline(x=1/per, color='r', linestyle='dashed')        
            ax9.plot(freqs, power, '.k', ms=1)
            ax9.set_xlabel('Frequency [1/days]')
            ax9.set_xlim(ax7.get_xlim())

            if bls:
                ax6.set_ylabel('ZTF BLS Power')        
                ax9.set_ylabel('ZTF BLS Power')
            else:
                ax6.set_ylabel('ZTF LS Power') 
                ax9.set_ylabel('ZTF LS Power')   
            
        folded_t, folded_y, folded_dy = bin_timeseries(t%per_ztf, y, bins, dy=dy)
        plot_phase_curve(ax2, folded_t, folded_y, folded_dy, period=per_ztf,
                         ylabel="ZTF Relative Flux")


    # -- hr diagram ------------------------------------------------------------

    hr_diagram(gaia_tab, ra, dec, ax3)
    
    fig.tight_layout()
    plt.savefig(out_dir+suffix+'_panel.png', dpi=300)
    print('Saved '+out_dir+suffix+'_panel.png')
    plt.close()

