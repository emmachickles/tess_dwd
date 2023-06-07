# -- inputs --------------------------------------------------------------------

n_std = 5
detrend = "wotan"
wind = 0.1
pmin = 400 / 60 
# pmax = 0.15
pmax = 10
qmin = 0.01
qmax = 0.15

# wid_threshold=6
# pow_threshold=25
wid_threshold = 0
pow_threshold = 0
objid_type = None

wd_tab= "/scratch/echickle/WDs.txt"
wd_main = "/scratch/echickle/GaiaEDR3_WD_main.fits"
rp_ext = "/scratch/echickle/GaiaEDR3_WD_RPM_ext.fits"
qflag_dir = "/scratch/echickle/QLPqflags/"

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
    ticid = np.int64(ticid)
    suffix = "-{}-{}.npy".format(cam, ccd)
    # y = np.load(data_dir+'lc'+suffix)
    y = np.load(data_dir+'lc-2-2-ap1.1-in1.8-out2.3.npy')
    t = np.load(data_dir+'ts'+suffix)
    coord = np.load(data_dir+'co'+suffix)
    ticid_ccd = np.load(data_dir+'id'+suffix).astype('int')
    ind = np.nonzero(ticid_ccd == ticid)[0][0]
    y, coord = y[ind], coord[ind]
    ra, dec = coord[0], coord[1]
    print('Loaded S{}-{}-{} TIC{}'.format(sector, cam, ccd, ticid))

    # >> remove nonzero quality flags
    sector_dir = qflag_dir + 'sec%d/' % sector
    cn = np.load(data_dir+'cn'+suffix)
    file_names = os.listdir(sector_dir)
    file_names = [f for f in file_names if 'cam%dccd%d'%(cam, ccd) in f]
    qflag_data = []
    for f in file_names:
        qflag_data.extend(np.loadtxt(sector_dir+f))
    qflag_data = np.array(qflag_data)
    bad_inds = np.nonzero(qflag_data[:,1])[0]
    bad_cadence = qflag_data[:,0][bad_inds]
    _, comm1, comm2 = np.intersect1d(cn, bad_cadence, return_indices=True)

    cn = np.delete(cn, comm1)
    t = np.delete(t, comm1)
    y = np.delete(y, comm1)

    # >> detrend and remove outliers
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
                           objid=ticid, objid_type=objid_type, suffix=suffix, 
                           ra=ra, dec=dec,
                           pow_threshold=pow_threshold, wid_threshold=wid_threshold,
                           wd_main=wd_main, rp_ext=rp_ext, wd_tab=wd_tab)
        print(res)
        return [ticid, ra, dec] + list(res)

    except:
        print('!! Failed !!')
        res = [ticid, ra, dec] + list(np.zeros(12))
        return res

if __name__ == '__main__':

    import sys, pdb
    sector, cam, ccd = 61, 1, 1
    ticid = 178366477
    data_dir = "/scratch/data/tess/lcur/ffi/s0061-lc-ZTF/"
    output_dir = "/scratch/echickle/s%04d-ZTF/"%sector
    bls_dir = output_dir + "s%04d"%sector + "-bls-{}-{}/".format(cam,ccd)

    p = [sector, cam, ccd, ticid, data_dir, bls_dir]
    res = run_process(p)
    print(res)
