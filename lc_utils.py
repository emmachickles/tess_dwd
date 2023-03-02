import numpy as np
import matplotlib.pyplot as plt
import pdb

def extract_lc():
    import pandas as pd
    
    # >> load white dwarf catalog
    wd_cat  = pd.read_csv('/home/echickle/work/tess_dwd/WDs.txt', header=None, sep='\s+')

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

def bin_timeseries(t, y, dy, bins):
    trunc = len(y) // bins * bins
    t_binned = np.split(t[:trunc], bins)
    y_binned = np.split(y[:trunc], bins)
    dy_binned = np.split(dy[:trunc], bins)

    if trunc < len(y):
        y_binned[-1] = np.append(y_binned[-1], y[trunc:])
        t_binned[-1] = np.append(t_binned[-1], t[trunc:])
        dy_binned[-1] = np.append(dy_binned[-1], dy[trunc:])

    for i in range(bins):
        t_binned[i] = np.average(t_binned[i])

        # dy_binned[i] = 1 / np.sqrt(np.sum(dy_binned[i]))
        err = np.std(y_binned[i])
        
        y_binned[i] = np.average(y_binned[i], weights=dy_binned[i])
        dy_binned[i] = err

    return t_binned, y_binned, dy_binned


def prep_lc(t, y, n_std=2, wind=0.1, lim=1000, diag=False, ticid=None, cam=None,
            ccd=None, coord=None, output_dir=None):
    from wotan import flatten

    # if diag:
    #     import sys
    #     sys.path.insert(0, "/home/submit/echickle/work/")    
    #     from KBB_Utils.KBB_Utils.Period_Finding import BLS
    from Period_Finding import BLS

    flag = False
    
    # >> remove nans
    inds = np.nonzero(~np.isnan(y))
    t, y = t[inds], y[inds]

    # >> sigma-clip         
    med = np.median(y)
    std = np.std(y)
    inds = np.nonzero( (y > med - n_std*std) * (y < med + n_std*std) )
    t, y = t[inds], y[inds]

    if diag:
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
        
    
    # >> normalize        
    if np.median(y) < 0:
        y = y / np.abs(med) + 2.
    else:
        y = y / med
    
    # >> detrending 
    y = flatten(t, y, window_length=wind, method='biweight')
    inds = np.nonzero(~np.isnan(y))
    t, y = t[inds], y[inds]

    if diag:
        dy = np.ones(y.shape)                
        t, y, dy, period, bls_power_best, freqs, power, dur, epo = \
            BLS(t,y,dy,pmin=7,pmax=0.25,qmin=0.005,qmax=0.2,remove=True)
        prefix = 'TIC%016d'%ticid+'_1_cam_'+str(cam)+'_ccd_'+str(ccd)+\
            '_pow_'+str(bls_power_best)+'_per_'+str(round(period*1440,2))+\
            '_dur_'+str(dur)+'_epo_'+str(epo)+\
            '_ra_{}_dec_{}_'.format(coord[0], coord[1])                
        make_phase_curve(t, y, period, dy=dy, output_dir=output_dir,
                             prefix=prefix, freqs=freqs, power=power,
                             ticid=ticid, bins=100)
    
    if y.shape[0] < lim:
        flag = True
    
    
    # if diag:
    #     make_phase_curve(t, y, 49.71/1440, dy=np.ones(y.shape)*0.1,
    #                      output_dir='/home/submit/echickle/foo/',
    #                      prefix='nstd'+str(n_std)+'-wind'+str(wind)+'-2-')
        
    return t, y, flag
    

def make_phase_curve(t, y, period, dy=None, output_dir=None, prefix='', freqs=None, power=None, ticid=None,
                     bins=200, t0=None, bls=True):
    '''TODO:
    * bin phase curve
    * change x axis to orbital phase
    * repeat phase curve twice '''
    from astropy.timeseries import TimeSeries
    from astropy.time import Time
    import astropy.units as u

    # >> normalize
    med = np.median(y)
    if np.median(y) < 0:
        y = y / np.abs(med) + 2.
    else:
        y = y/med    
    
    if type(dy) == type(None):
        dy = np.ones(y.shape)

    # >> fold
    t = t % period 
    inds = np.argsort(t)
    folded_t, folded_y, folded_dy = t[inds], y[inds], dy[inds]
        
    folded_t, folded_y, folded_dy = bin_timeseries(folded_t, folded_y, folded_dy, bins)
    if output_dir:
        if type(freqs) == type(None):
            _,_,_, _, bls_power_best, freqs, power, dur, epo = \
                BLS(t,y,dy,pmin=7,pmax=0.25,qmin=0.005,qmax=0.2,remove=False)

            
        fig, ax = plt.subplots(nrows=3, figsize=(8, 8))
        if type(ticid) != type(None):
            ax[0].set_title('TIC '+str(ticid)+'\nperiod: '+str(round(period*1440,2))+' min')
        else:
            ax[0].set_title('period: '+str(round(period*1440,2))+' min')
        ax[0].plot(freqs, power, '.k', ms=1, alpha=0.5)
        ax[0].set_xlabel('Frequency [1/days]')

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
        if bins:
            ax[2].errorbar(folded_t*1440, folded_y, yerr=folded_dy*0.1,
                           fmt='.k', ms=1, elinewidth=1)
            ax[2].errorbar((folded_t+shift)*1440, folded_y, yerr=folded_dy*0.1,
                           fmt='.k', ms=1, elinewidth=1)

        else:
            ax[2].plot(folded_t*1440, folded_y, '.k', ms=1)
            ax[2].plot((folded_t+shift)*1440, folded_y, '.k', ms=1)            
        ax[2].set_xlabel('Time [minutes]')
        ax[2].set_ylabel('Relative Flux')
        fig.tight_layout()

        plt.savefig(output_dir+prefix+'phase_curve.png', dpi=300)
        print('Saved '+output_dir+prefix+'phase_curve.png')
        plt.close()

    else:
        return np.array(folded_t), np.array(folded_y), np.array(folded_dy)

def plot_phase_curve(t, y, per, out_dir, bins=200):
    inds = np.nonzero(~np.isnan(y))
    t, y = t[inds], y[inds]
    folded_t, folded_y, folded_dy = make_phase_curve(t, y, per, bins=bins)
    shift = np.max(folded_t) - np.min(folded_t)
    fig, ax = plt.subplots(figsize=(10,6))
    if type(bins) != type(None):
        ax.errorbar(folded_t*1440, folded_y, yerr=folded_dy*0.1,
                       fmt='.k', ms=1, elinewidth=1)
        ax.errorbar((folded_t+shift)*1440, folded_y, yerr=folded_dy*0.1,
                       fmt='.k', ms=1, elinewidth=1)
    else:
        ax.plot(folded_t*1440, folded_y, '.k', ms=1, alpha=0.3)
        ax.plot((folded_t+shift)*1440, folded_y, '.k', ms=1, alpha=0.3) 
    ax.set_xlabel('Time [minutes]')
    ax.set_ylabel('Relative Flux')
    fig.tight_layout()
    plt.savefig(out_dir+'phase_curve.png', dpi=300)
    print('Saved '+out_dir+'phase_curve.png')

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

def make_panel_plot(t,y,freqs,power,period,prefix, bins=200):
    gs = fig.add_gridspec(nrows=2, ncols=2, figsize=(10,10))
    ax1 = fig.add_subplot(gs[:, 0])
    ax2 = fig.add_subplot(gs[0, 1])
    ax3 = fig.add_subplot(gs[1, 1])
    
    inds = np.nonzero(~np.isnan(y))
    t, y = t[inds], y[inds]
    folded_t, folded_y, folded_dy = make_phase_curve(t, y, per, bins=bins)
    shift = np.max(folded_t) - np.min(folded_t)
    fig, ax = plt.subplots(figsize=(10,6))
    if type(bins) != type(None):
        ax.errorbar(folded_t*1440, folded_y, yerr=folded_dy*0.1,
                       fmt='.k', ms=1, elinewidth=1)
        ax.errorbar((folded_t+shift)*1440, folded_y, yerr=folded_dy*0.1,
                       fmt='.k', ms=1, elinewidth=1)
    else:
        ax.plot(folded_t*1440, folded_y, '.k', ms=1, alpha=0.3)
        ax.plot((folded_t+shift)*1440, folded_y, '.k', ms=1, alpha=0.3) 
    ax.set_xlabel('Time [minutes]')
    ax.set_ylabel('Relative Flux')
    fig.tight_layout()
    plt.savefig(out_dir+'phase_curve.png', dpi=300)
    print('Saved '+out_dir+'phase_curve.png')
    plt.close()
