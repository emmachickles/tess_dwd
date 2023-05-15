import lc_utils as lcu
from Period_Finding import BLS
import os
import pdb
import matplotlib.pyplot as plt
import numpy as np

pos_iqr = 3
neg_iqr = 20
skiprows = 1
objid_type = None

qmin = 0.01
qmax = 0.15

bins=100

output_dir = "/home/echickle/out/DWD_pc/"
os.makedirs(output_dir, exist_ok=True)
data_dir = "/home/echickle/data/atlasforcedphotometryresults_DWD/"
wd_main = "/home/echickle/data/GaiaEDR3_WD_main.fits"
rp_ext = "/home/echickle/data/GaiaEDR3_WD_RPM_ext.fits"

# period = [0.0800126, 0.09986, 0.11601549, 0.2350606, 0.246137]
period = [0.11601549, 0.0800126, 0.2350606, 0.09986, 0.246137]

fnames = os.listdir(data_dir)
for i in range(len(fnames)):
    f = fnames[i]

    t, y, dy, ra, dec = lcu.load_atlas_lc(data_dir+f, pos_iqr=pos_iqr, neg_iqr=neg_iqr,
                                          skiprows=skiprows, clip=False)

    print(ra)
    print(dec)
    gaiaid = f.split('/')[-1]
    suffix = '_GID_'+gaiaid+'_ra_{}_dec_{}_'.format(ra, dec)    

    # t, y, dy, period, bls_power_best, freqs, power, q, phi0 = \
    #     BLS(t,y,dy,pmin=2,pmax=10,qmin=qmin,qmax=qmax,remove=False)
    # res = lcu.vet_plot(t, y, freqs, power, q, phi0, output_dir=output_dir,
    #                    objid=gaiaid, objid_type=objid_type, ra=ra, dec=dec,
    #              dy=dy, suffix=suffix, wd_main=wd_main, rp_ext=rp_ext)
    # # per, q, epo = res[3], res[5], res[7]
    # lcu.plot_eclipse_timing(t, y, per, epo, q, output_dir+'GAIAID_{}_{}_{}_'.format(gid, ra, dec))

    per = period[i]
    y, dy = lcu.normalize_lc(y, dy)

    fig, ax = plt.subplots()
    folded_t, folded_y, folded_dy = lcu.bin_timeseries(t%per, y, bins, dy=dy)
    lcu.plot_phase_curve(ax, folded_t, folded_y, folded_dy, period=per)
    
    # w = int(0.1*len(y))
    # tconv = np.convolve(t%per, np.ones(w), 'valid') / w
    # yconv = np.convolve(y, np.ones(w), 'valid') / w
    # inds = np.argsort(tconv)
    # tconv, yconv = tconv[inds], yconv[inds]    
    # ax.plot(tconv*1440, yconv, '-')

    # shift = np.max(t%per) - np.min(t%per)
    # ax.plot(t%per*1440, y, '.k', ms=1)
    # ax.plot((t%per + shift)*1440, y, '.k', ms=1)    
    # ax.text(0.95, 0.05, str(np.round(per*1440,5))+" min",
    #         horizontalalignment="right", transform=ax.transAxes)
    # ax.set_xlabel('Time [minutes]')
    # ax.set_ylabel('Relative flux')
    
    fig.savefig(output_dir + 'ra_{}_dec_{}_phase_curve.png'.format(ra, dec))
    
