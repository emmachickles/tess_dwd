import pdb
from vet_utils import *
from lc_utils import get_tess_lc, normalize_lc, bin_timeseries, plot_phase_curve, prep_lc
import numpy as np
import matplotlib.pyplot as plt
import os

out_dir = '/scratch/echickle/cut/'
os.makedirs(out_dir, exist_ok=True)
data_dir = '/scratch/echickle/tess/BLS_results/'
tess_dir = "/scratch/data/tess/lcur/ffi/"
sector_list = [56,57,58,59,60,61,62,63,64,65]
bins=100
figsize=(5,2)

dur_min = 3 # min
dur_max = 8 # min
per_min = 60
per_max = 120 # min
nt_min = 50
dphi_max = 0.03
pow_min = 25
wid_min = 6
dec_max = -25

# Load Gaia white dwarf results 
# ticid, ra, dec, power, snr, wid, per, nt, dphi, dur
result_list = append_result_file(data_dir, sector_list)

# Apply cut
good_idx = np.nonzero( (result_list[:,6]<per_max/1440.) * \
                       (result_list[:,6]>per_min/1440.) *\
                       (result_list[:,7]>nt_min) * \
                       (result_list[:,8]<dphi_max) * \
                       (result_list[:,3]>pow_min) * \
                       (result_list[:,5]>wid_min) * \
                       (result_list[:,9]>dur_min) *\
                       (result_list[:,9]<dur_max) * \
                       (result_list[:,2]<dec_max))

ra_list = result_list[:,1][good_idx]
dec_list = result_list[:,2][good_idx]
per_list=result_list[:,6][good_idx]
for i in range(len(good_idx[0])):
    # ra, dec = 121.174886, -2.262535
    ra = result_list[:,1][good_idx][i]
    dec = result_list[:,2][good_idx][i]
    per = result_list[:,6][good_idx][i]
    dur = result_list[:,9][good_idx][i]
    nt = result_list[:,7][good_idx][i]
    dphi = result_list[:,8][good_idx][i]
    wid = result_list[:,5][good_idx][i]

    suffix='ra_{}_dec_{}_per_{}_dur_{}_nt_{}_dphi_{}_wid_{}'.format(ra,dec, per,dur,nt,dphi,wid)
    ticid, sector, cam, ccd = get_tess_lc(tess_dir, ra=ra, dec=dec, uzay=True)
    ccd_dir = tess_dir+'s00{}-lc/'.format(sector)
        
    f_suffix = '-{}-{}.npy'.format(cam, ccd)        
    t = np.load(ccd_dir+'ts'+f_suffix)
    ticid_list = np.load(ccd_dir+'id'+f_suffix)
    ticid, ticid_list = np.int64(ticid), np.int64(ticid_list)
    ind = np.nonzero(ticid_list == ticid)[0][0]
    y = np.load(ccd_dir+'lc'+f_suffix)[ind]
    t, y, flag = prep_lc(t, y)
    
    y, _ = normalize_lc(y)
    folded_t, folded_y, folded_dy = bin_timeseries(t%per, y, bins)

    
    fig00, ax00=plt.subplots(figsize=figsize)
    folded_t, folded_y, folded_dy = bin_timeseries(t%per, y, bins)
    plot_phase_curve(ax00, folded_t, folded_y, folded_dy, period=per,
                     ylabel="TESS Relative Flux")
    fig00.savefig(out_dir+suffix+'_binned_TESS.png', dpi=300)
    print('Saved '+out_dir+suffix+'_binned_TESS.png')
    plt.close()

pdb.set_trace()
