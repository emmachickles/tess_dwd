# -- inputs -------------------------------------------------------------------

pmin = 2 # minutes
pmax = 10 # days 
qmin = 0.01
qmax = 0.15
dlogq = 0.1

pow_threshold=25
snr_threshold=1.
per_threshold=210
wid_threshold=5

# data_dir = "/home/echickle/ATLAS_TEST/"
# data_dir = "/matchfiles/data2/ATLAS/"
# output_dir = "/home/echickle/out/"
wd_main = "/data/GaiaEDR3_WD_main.fits"
rp_ext = "/data/GaiaEDR3_WD_RPM_ext.fits"

data_dir = "/pool001/echickle/ATLAS/"
output_dir = "/pool001/echickle/out/"
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
                       pow_threshold=pow_threshold, per_threshold=per_threshold,
                       snr_threshold=snr_threshold, wid_threshold=wid_threshold,
                       wd_main=wd_main, rp_ext=rp_ext)
    # res: sig, snr, wid, period, period_min, q, phi0, dur, epo, epo, rp, nt, dphi
    # end=time.time()
    print('Postprocessed light curves')
    # print(end-start)
    res = [int(gaiaid)] + res
    
    return res
    

if __name__ == '__main__':

    import sys
    fname = sys.argv[1]

    import os
    os.makedirs(output_dir, exist_ok=True)
    os.makedirs(bls_dir, exist_ok=True)

    # --------------------------------------------------------------------------

    run_process(data_dir + fname)


    
