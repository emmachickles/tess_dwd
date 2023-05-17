# -- inputs --------------------------------------------------------------------

n_std = 5
detrend = "wotan"
wind = 0.1
pmin = 400 / 60 
# pmax = 0.15
qmin = 0.01
qmax = 0.15

wid_threshold=6
pow_threshold=25

def run_process(p):
    sector, cam, ccd, ticid, data_dir, bls_dir = p

    import numpy as np
    import lc_utils as lcu
    from Period_Finding import BLS
    import os
    import matplotlib as mpl
    mpl.use('Agg')

    print('Starting S{}-{}-{} TIC{}'.format(sector, cam, ccd, ticid))
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
        res = [ticid, coord[0], coord[1]] + list(np.zeros(12))
        return res

    try:
        dy = np.ones(y.shape)
        freqs_to_remove = []

        df = 0.1
        freqs_to_remove.append([86400/(200*2) - df, 86400/(200*2) + df])
        freqs_to_remove.append([86400/500 - df, 86400/500 + df])    
        freqs_to_remove.append([86400/(200*3) - df, 86400/(200*3) + df])
        freqs_to_remove.append([86400/600 - df, 86400/600 + df])    
        freqs_to_remove.append([86400/(200*4) - df, 86400/(200*4) + df])
        freqs_to_remove.append([86400/(200*5) - df, 86400/(200*5) + df])     
        freqs_to_remove.append([86400/(200*6) - df, 86400/(200*6) + df]) 
        freqs_to_remove.append([86400/(200*7) - df, 86400/(200*7) + df])   

        
        t, y, dy, period, bls_power_best, freqs, power, q, phi0 = \
            BLS(t,y,dy,pmin=pmin,pmax=pmax,qmin=qmin,qmax=qmax,freqs_to_remove=freqs_to_remove)
        print('Computed BLS S{}-{}-{} TIC{}'.format(sector, cam, ccd, ticid))
        
        suffix = '_TIC%016d'%ticid+'_s%04d_'%sector+'cam_'+\
                 str(cam)+'_ccd_'+str(ccd)+\
                 '_ra_{}_dec_{}_'.format(coord[0], coord[1])
        res = lcu.vet_plot(t, y, freqs, power, q, phi0, output_dir=bls_dir,
                           objid=ticid, 
                           suffix=suffix, bins=100, save_npy=False,
                           pow_threshold=pow_threshold, wid_threshold=wid_threshold)
        return [ticid, coord[0], coord[1]] + list(res)

    except:
        res = [ticid, coord[0], coord[1]] + list(np.zeros(12))
        return res

if __name__ == '__main__':

    import sys
    sector, cam, ccd = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
    ticid = int(sys.argv[4])

    p = [sector, cam, ccd, ticid]
    run_process(p)
