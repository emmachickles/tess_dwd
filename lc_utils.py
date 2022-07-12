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
    
def make_phase_curve(t, y, period, plot=True, output_dir='', prefix='', freqs=None, power=None, ticid=None):
    '''TODO:
    * bin phase curve
    * change x axis to orbital phase
    * repeat phase curve twice '''
    from astropy.timeseries import TimeSeries
    from astropy.time import Time
    import astropy.units as u

    ts = TimeSeries(time=Time(t, format='jd'), data={'flux': y})
    ts_folded = ts.fold(period=period*u.d) 

    sorted_inds = np.argsort(ts_folded['time'])
    folded_t = ts_folded['time'][sorted_inds]
    folded_y = ts_folded['flux'][sorted_inds]

    if plot:
        fig, ax = plt.subplots(nrows=2)
        ax[0].set_title('TIC '+str(ticid)+'\nperiod: '+str(round(period*1440,2))+' min')
        ax[0].plot(freqs, power, '.k', ms=1)
        ax[0].set_xlim([0,50])
        ax[0].set_xlabel('Frequency [1/days]')
        ax[0].set_ylabel('BLS Power')
        ax[1].plot(folded_t.value, folded_y.value, '.k', ms=1)
        ax[1].set_xlabel('Phase')
        ax[1].set_ylabel('PDCSAP_FLUX')
        # fname = '/data/submit/echickle/out/TIC1201247611_pc_qmin'+str(qmin)+'_qmax'+str(qmax)+'_period'+str(period)+'_sig'+str(bls_power_best)+'.png'
        # fname = '/data/submit/echickle/out/TIC1201247611-%02d'%n_iter+'.png'
        # fname = img_dir+'TIC1201247611-%02d'%n_iter+'-p'+str(period)+'.png'
        fig.tight_layout()



        # plt.figure()
        # norm_y = folded_y.value / np.nanmedian(folded_y.value)
        # plt.plot(folded_t.jd, norm_y, '.k', ms=1)
        # plt.xlabel('Time [days]')
        # plt.ylabel('Flux')
        plt.savefig(output_dir+prefix+'phase_curve.png')
        print('Saved '+output_dir+prefix+'phase_curve.png')

    return
