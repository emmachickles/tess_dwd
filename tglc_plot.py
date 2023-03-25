import os
from astropy.io import fits
import lc_utils as lcu
import pdb
import numpy as np
import matplotlib.pyplot as plt
import sys
from Period_Finding import BLS, LS_Astropy

n_std=5
detrend = "polyfit"
pmin = 410 / 60
pmax = 0.13
qmin = 0.01
qmax = 0.15

# hdul = fits.open('/data/submit/echickle/foo/lc/hlsp_tglc_tess_ffi_gaiaid-1629388752470472704-s0014-cam3-ccd1_tess_v1_llc.fits')
hdul = fits.open('/data/submit/echickle/foo1/lc/hlsp_tglc_tess_ffi_gaiaid-2239471475135041664-s0060-cam4-ccd4_tess_v1_llc.fits')
# hdul = fits.open('/data/submit/echickle/foo1/lc/hlsp_tglc_tess_ffi_gaiaid-1629388752470472704-s0041-cam3-ccd4_tess_v1_llc.fits')

t = hdul[1].data['time']
y = hdul[1].data['cal_psf_flux']
per = 25./1440
ticid = 1884076062

output_dir = "/data/submit/echickle/amcvn_test/"
os.makedirs(output_dir, exist_ok=True)
bls_dir    = output_dir+"bls/"
ls_dir     = output_dir+"ls/"
os.makedirs(bls_dir, exist_ok=True)
os.makedirs(ls_dir, exist_ok=True)


# per = 49.7/1440
# per = 1/20.727432959364357

t, y, flag = lcu.prep_lc(t, y, n_std=n_std, detrend=detrend)
dy = np.ones(y.shape) * 0.1
# # >> fold
# t = t % per * 1440
# inds = np.argsort(t)
# t, y = t[inds], y[inds]

# # >> bin
# dy = np.ones(y.shape)
# t, y, dy = lcu.bin_timeseries(t, y, 60, dy=dy)

# # >> plot
# shift = np.max(t) - np.min(t)
# plt.figure(figsize=(5,3))
# plt.errorbar(t, y, yerr=np.ones(len(y))*0.01,
#                fmt='.k', ms=1, elinewidth=1)
# plt.errorbar(t+shift, y, yerr=np.ones(len(y))*0.01,
#                fmt='.k', ms=1, elinewidth=1)
# plt.xlabel('Time [minutes]')
# plt.ylabel('Relative Flux')
# plt.tight_layout()
# plt.savefig(output_dir+'TIC'+str(ticid)+'_tglc_phase_curve.png', dpi=300)
# print(output_dir+'TIC'+str(ticid)+'_tglc_phase_curve.png')

_, _, _, per, pow_best, bls_frq, bls_pow, dur, epo, delta = \
        BLS(t,y,dy,pmin=pmin,pmax=pmax,qmin=qmin,qmax=qmax,remove=True)
_, _, _, ls_per, ls_pow_best, ls_frq, ls_pow = LS_Astropy(t,y,dy,pmax=pmax)

prefix = ni+'pow_'+str(pow_best)+'_delta_'+str(round(delta,5))+\
         '_per_'+str(round(per*1440,5))+\
    '_TIC%016d'%ticid+'_dur_'+str(dur)+'_epo_'+str(epo)+\
    '_ra_{}_dec_{}_'.format(coord[0], coord[1])                

lcu.make_phase_curve(t, y, per, output_dir=bls_dir, prefix=prefix,
                     freqs=bls_frq, power=bls_pow, ticid=ticid, bins=100,
                     dur=dur, epo=epo)

prefix = ni+'pow_'+str(ls_pow_best)+'_per_'+str(round(ls_per*1440,5))+\
    '_TIC%016d'%ticid+'_ra_{}_dec_{}_'.format(coord[0], coord[1])                

lcu.make_phase_curve(t, y, ls_per, output_dir=ls_dir, prefix=prefix,
                     freqs=ls_frq, power=ls_pow, ticid=ticid, bins=100,
                     bls=False)


