# -- inputs -------------------------------------------------------------------

pos_iqr = 3
neg_iqr = 10
# skiprows = 1
objid_type = None
skiprows = 0
# objid_type = 'GAIAID'

pmin = 2 # minutes
pmax = 10 # days 
qmin = 0.01
qmax = 0.15
dlogq = 0.1

pow_threshold=10.
# pow_threshold=25
snr_threshold=4.
per_threshold=210
wid_threshold=5

# wd_main = "/home/echickle/data/GaiaEDR3_WD_main.fits"
# rp_ext = "/home/echickle/data/GaiaEDR3_WD_RPM_ext.fits"
# wd_main = "/data/GaiaEDR3_WD_main.fits"
# rp_ext = "/data/GaiaEDR3_WD_RPM_ext.fits"
wd_main = "/nobackup1c/users/echickle/GaiaEDR3_WD_main.fits"
rp_ext = "/nobackup1c/users/echickle/GaiaEDR3_WD_RPM_ext.fits"

def run_process(p):
    import lc_utils as lcu
    from Period_Finding import BLS
    import matplotlib as mpl
    mpl.use('Agg')
    # import time

    f, data_dir, bls_dir = p

    print('Starting '+f)
    
    t, y, dy, ra, dec = lcu.load_atlas_lc(f, pos_iqr=pos_iqr, neg_iqr=neg_iqr,
                                          skiprows=skiprows)
    print('Loaded '+f)

    # start =time.time()
    freqs_to_remove = []
    df = 0.05
    freqs_to_remove.append([1 - df, 1 + df])
    freqs_to_remove.append([1/2. - df, 1/2. + df])
    freqs_to_remove.append([1/4. - df, 1/4. + df])    
    y, dy = lcu.normalize_lc(y, dy)
    t, y, dy, period, bls_power_best, freqs, power, q, phi0 = \
        BLS(t,y,dy,pmin=pmin,pmax=pmax,qmin=qmin,qmax=qmax,dlogq=dlogq,
            freqs_to_remove=freqs_to_remove)
    # end=time.time()
    print('Computed BLS '+f)        
    # print(end-start)
    gaiaid = f.split('/')[-1]
    suffix = '_GID_'+gaiaid+'_ra_{}_dec_{}_'.format(ra, dec)    

    # start = time.time()
    res = lcu.vet_plot(t, y, freqs, power, q, phi0, output_dir=bls_dir, dy=dy,
                       suffix=suffix, objid=gaiaid, objid_type=objid_type,
                       pow_threshold=pow_threshold, per_threshold=per_threshold,
                       snr_threshold=snr_threshold, wid_threshold=wid_threshold,
                       wd_main=wd_main, rp_ext=rp_ext, ra=ra, dec=dec)
    # res: sig, snr, wid, period, period_min, q, phi0, epo, rp, nt, dphi
    # end=time.time()
    print('Postprocessed light curves')
    # print(end-start)
    res = [gaiaid, ra, dec] + list(res)
    
    return res
    

if __name__ == '__main__':


    import os
    data_dir = "/home/echickle/data/atlasforcedphotometryresults_DWD/"
    output_dir = "/home/echickle/out/"
    bls_dir = output_dir + 'DWD_plot/'    
    os.makedirs(output_dir, exist_ok=True)
    os.makedirs(bls_dir, exist_ok=True)

    import sys
    # fname = sys.argv[1]
    fname = data_dir+'job474720.txt'
    
    # --------------------------------------------------------------------------

    p = [fname, data_dir, bls_dir]
    
    run_process(p)


    
