import matplotlib.pyplot as plt
import os
import pdb
import numpy as np
from Period_Finding import BLS
import lc_utils as lcu

data_dir = '/home/echickle/data/atlasforcedphotometryresults_DWD/'
fnames = sorted(os.listdir(data_dir))
output_dir = '/home/echickle/out/combine_test/'
os.makedirs(output_dir, exist_ok=True)
per_true = [0.0800126, 0.09986, 0.11601549, 0.2350606, 0.246137]
home_dir = '/home/echickle/'

clip = False
pos_iqr = 3
neg_iqr = 10
skiprows = 1
objid_type = None
pmin = 2 # minutes
pmax = 10 # days 
qmin = 0.01
qmax = 0.15
dlogq = 0.1
freqs_to_remove = []
df = 0.05
freqs_to_remove.append([1 - df, 1 + df])
freqs_to_remove.append([1/2. - df, 1/2. + df])
freqs_to_remove.append([1/4. - df, 1/4. + df])    

for i in [2,4]:
    fname_atlas = data_dir+fnames[i]
    t_atlas, y_atlas, dy_atlas, ra, dec  = lcu.load_atlas_lc(fname_atlas, pos_iqr=pos_iqr,
                                                             neg_iqr=neg_iqr, skiprows=skiprows,
                                                             clip=clip)
    print('RA {} Dec {}'.format(ra, dec))
    t_atlas, y_atlas, dy_atlas, period, bls_power_best_atlas, freqs, power, q, phi0 = \
        BLS(t_atlas,y_atlas,dy_atlas,pmin=pmin,pmax=pmax,qmin=qmin,qmax=qmax,dlogq=dlogq,
            freqs_to_remove=freqs_to_remove)
    suffix = '_ATLAS_ra_{}_dec_{}_per_{}'.format(ra, dec, round(per_true[i]*1440,2))
    lcu.plot_signal(t_atlas, y_atlas, dy_atlas, freqs, power, per_true[i], output_dir, suffix=suffix)


    # fnames_ztf = lcu.get_ztf_lc(ra, dec)
    fnames_ztf = []
    for filt in ['g', 'i', 'r']:
        f = home_dir + '{}_{}_{}.lc'.format(ra, dec, filt)
        if os.path.exists(f):
            if os.stat(f).st_size > 0:
                fnames_ztf.append(f)
    t_ztf, y_ztf, dy_ztf = lcu.load_ztf_lc(fnames_ztf)
    y_ztf, dy_ztf = lcu.normalize_lc(y_ztf, dy_ztf)
    t_ztf, y_ztf, dy_ztf, period, bls_power_best_ztf, freqs, power, q, phi0 = \
        BLS(t_ztf,y_ztf,dy_ztf,pmin=pmin,pmax=pmax,qmin=qmin,qmax=qmax,dlogq=dlogq,
            freqs_to_remove=freqs_to_remove)
    suffix = '_ZTF_ra_{}_dec_{}_per_{}'.format(ra, dec, round(per_true[i]*1440,2))
    lcu.plot_signal(t_ztf, y_ztf, dy_ztf, freqs, power, per_true[i], output_dir, suffix=suffix)
    

    # Combine ATLAS and ZTF light curves
    y_ztf += np.median(y_atlas) - np.median(y_ztf)
    
    t_comb = np.append(t_atlas, t_ztf)
    y_comb = np.append(y_atlas, y_ztf)
    dy_comb = np.append(dy_atlas, dy_ztf)
    inds = np.argsort(t_comb)
    t_comb, y_comb, dy_comb=  t_comb[inds], y_comb[inds], dy_comb[inds]
    t_comb, y_comb, dy_comb, period, bls_power_best_ztf, freqs, power, q, phi0 = \
        BLS(t_comb,y_comb,dy_comb,pmin=pmin,pmax=pmax,qmin=qmin,qmax=qmax,dlogq=dlogq,
            freqs_to_remove=freqs_to_remove)
    suffix = '_COMBINED_ra_{}_dec_{}_per_{}'.format(ra, dec, round(per_true[i]*1440,2))
    lcu.plot_signal(t_comb, y_comb, dy_comb, freqs, power, per_true[i], output_dir, suffix=suffix)
