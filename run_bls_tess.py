# -- inputs --------------------------------------------------------------------

n_std = 5
detrend = "wotan"
wind = 0.1
pmin = 400 / 60 
# pmax = 0.15
pmax = 13
qmin = 0.01
qmax = 0.15

def run_process(p):
    sector, cam, ccd, ticid = p

    import numpy as np
    import lc_utils as lcu
    from Period_Finding import BLS
    import os
    import matplotlib as mpl
    mpl.use('Agg')

    print('Starting S{}-{}-{} TIC{}'.format(sector, cam, ccd, ticid))

    data_dir = "/scratch/data/tess/lcur/ffi/s%04d-lc/"%sector
    # data_dir = "/scratch/data/tess/lcur/ffi/s%04d-lc-ZTF/"%sector

    output_dir = "/scratch/echickle/s%04d/"%sector
    # output_dir = "/scratch/echickle/s%04d-ZTF/"%sector

    bls_dir = output_dir + "s%04d"%sector + "-bls-{}-{}/".format(cam,ccd)
    os.makedirs(output_dir, exist_ok=True)
    os.makedirs(bls_dir, exist_ok=True)

    # >> load data
    suffix = "-{}-{}.npy".format(cam, ccd)
    y = np.load(data_dir+'lc'+suffix)
    t = np.load(data_dir+'ts'+suffix)
    coord = np.load(data_dir+'co'+suffix)
    ticid_ccd = np.load(data_dir+'id'+suffix).astype('int')
    ind = np.nonzero(ticid_ccd == ticid)[0][0]
    y, coord = y[ind], coord[ind]
    print('Loaded S{}-{}-{} TIC{}'.format(sector, cam, ccd, ticid))

    t, y, flag = lcu.prep_lc(t, y, n_std=n_std, detrend=detrend, wind=wind)
    if flag:
        res = [ticid] + list(np.zeros(9, dtype='int'))
        return res

    try:
        dy = np.ones(y.shape)
        t, y, dy, period, bls_power_best, freqs, power, q, phi0 = \
            BLS(t,y,dy,pmin=pmin,pmax=pmax,qmin=qmin,qmax=qmax,remove=True)
        print('Computed BLS S{}-{}-{} TIC{}'.format(sector, cam, ccd, ticid))
        
        suffix = '_TIC%016d'%ticid+'_s%04d_'%sector+'cam_'+\
                 str(cam)+'_ccd_'+str(ccd)+\
                 '_ra_{}_dec_{}_'.format(coord[0], coord[1])
        res = lcu.vet_plot(t, y, freqs, power, q, phi0, output_dir=bls_dir,
                           objid=ticid, 
                           suffix=suffix, bins=100, save_npy=False)
        return [ticid] + list(res)

    except:
        res = [ticid] + list(np.zeros(9, dtype='int'))
        return res

if __name__ == '__main__':

    import sys
    sector, cam, ccd = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
    ticid = int(sys.argv[4])

    p = [sector, cam, ccd, ticid]
    run_process(p)
