def match_period(period, period_true, dt=2./1440):
    if (period > period_true-dt \
        and period < period_true+dt) or\
    (period > period_true*2-dt \
     and period < period_true*2+dt) or\
    (period > period_true/2-dt \
     and period < period_true/2+dt):
        return True
    else:
        return False

def read_result_file(fname, bls=True):
    import pandas as pd
    if bls:
        # ticid, ra, dec, sig, snr, wid, period, period_min, q, phi0, dur, epo, rp, nt, dphi
        cat = pd.read_csv(fname, header=None, skiprows=1)    
        ticid, ra, dec, power, snr, wid, per, nt, dphi = cat[0], cat[1], cat[2], cat[3], cat[4], cat[5], cat[6], cat[12], cat[13]
        ticid, wid, nt = np.int64(ticid), np.int64(wid), np.int64(nt)
        return ticid, ra, dec, power, snr, wid, per, nt, dphi
    else:
        # ticid, ra, dec, sig, wid, per, per_min, dphi
        cat = pd.read_csv(fname, header=None, skiprows=1)    
        ticid, ra, dec, power, wid, per, dphi = cat[0], cat[1], cat[2], cat[3], cat[4], cat[5], cat[7]
        ticid, wid = np.int64(ticid), np.int64(wid)
        return ticid, ra, dec, power, wid, per, dphi
        

def append_result_file(result_dir, sector_list, cam_list=[1,2,3,4], ccd_list=[1,2,3,4], bls=True):
    result_list = [np.empty()]
    for sector in sector_list:
        for cam in cam_list:
            for ccd in ccd_list:
                if bls:
                    fname = result_dir+'BLS-{}-{}-{}.result'.format(sector,cam,ccd)
                else:
                    fname = result_dir+'LS-{}-{}-{}.result'.format(sector,cam,ccd)
                result = read_result_file(fname, bls=bls)
                result = [np.array(arr) for arr in result]
                result_list.append(result)
    result_list = np.array(result_list)
    return
    
