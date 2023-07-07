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

def read_result_file(fname):
    # ticid, ra, dec, sig, snr, wid, period, period_min, q, phi0, dur, epo, rp, nt, dphi
    cat = pd.read_csv(fname, header=None, skiprows=1)    
    ticid, ra, dec, power, snr, wid, per, nt, dphi = cat[0], cat[1], cat[2], cat[3], cat[4], cat[5], cat[6], cat[12], cat[13]
    wid, nt = np.int64(wid), np.int64(nt)
    return ticid, ra, dec, power, snr, wid, per, nt, dphi

def append_result_file(mydir, sector_list, cam_list=[1,2,3,4], ccd_list=[1,2,3,4]):
    
    return
    
