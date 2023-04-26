# -- inputs -------------------------------------------------------------------

pmin = 2 # minutes
pmax = 10 # days 
qmin = 0.01
qmax = 0.15
dlogq = 0.1

# data_dir = "/home/echickle/ATLAS_TEST/"
data_dir = "/matchfiles/data2/ATLAS/"
output_dir = "/home/echickle/out/"

# data_dir = "/pool001/echickle/ATLAS/"
# output_dir = "/pool001/echickle/out/"
bls_dir    = output_dir + 'bls/'

def run_process(f):
    import lc_utils as lcu
    from Period_Finding import BLS
    import matplotlib as mpl
    mpl.use('Agg')
    # import time

    print('Starting '+f)
    
    t, y, dy, ra, dec = lcu.load_atlas_lc(f)
    print('Loaded '+f)

    # start =time.time()
    t, y, dy, period, bls_power_best, freqs, power, q, phi0 = \
        BLS(t,y,dy,pmin=pmin,pmax=pmax,qmin=qmin,qmax=qmax,dlogq=dlogq,remove=False)
    # end=time.time()
    print('Computed BLS '+f)        
    # print(end-start)
    gaiaid = f.split('/')[-1]
    suffix = '_GID_'+gaiaid+'_ra_{}_dec_{}_'.format(ra, dec)    

    # start = time.time()
    res = lcu.vet_plot(t, y, freqs, power, q, phi0, output_dir=bls_dir, dy=dy,
                       suffix=suffix, objid=gaiaid, objid_type='GAIAID',
                       plot_threshold=0) # !! 
    # end=time.time()
    print('Postprocessed light curves')
    # print(end-start)
    
    sig, snr, wid, period, period_min, q, phi0, dur, epo = res

    return int(gaiaid), sig, snr, wid, period, period_min, q, phi0, dur, epo
    

if __name__ == '__main__':

    import sys
    fname = sys.argv[1]

    import os
    os.makedirs(output_dir, exist_ok=True)
    os.makedirs(bls_dir, exist_ok=True)

    # --------------------------------------------------------------------------

    run_process(data_dir + fname)
