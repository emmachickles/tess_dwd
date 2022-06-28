import numpy as np
import matplotlib.pyplot as plt
import pdb

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
    
def make_phase_curve(t, y, period, plot=True, output_dir='', prefix=''):
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
        plt.figure()
        norm_y = folded_y.value / np.nanmedian(folded_y.value)
        plt.plot(folded_t.jd, norm_y, '.k', ms=1)
        plt.xlabel('Time [days]')
        plt.ylabel('Flux')
        plt.savefig(output_dir+prefix+'phase_curve.png')
        print('Saved '+output_dir+prefix+'phase_curve.png')

    return
