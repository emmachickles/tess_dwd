from lc_utils import get_tess_lc, normalize_lc, bin_timeseries, plot_phase_curve, prep_lc
import numpy as np
import matplotlib.pyplot as plt
import pdb
from Period_Finding import LS_Astropy
from astropy.timeseries import LombScargle

out_dir = "/home/echickle/out/vet/"
tess_dir = "/home/echickle/data/"
bins=100
ztf=False
figsize=(5,2)

ra_list = [231.90982999716]
dec_list = [-45.03539080131]
per_list = [745.877 / (1440*60.)]

for ra, dec, per in zip(ra_list, dec_list, per_list):
    # ra, dec = 121.174886, -2.262535
    suffix='ra_{}_dec_{}'.format(ra,dec)
    ticid, sector, cam, ccd = get_tess_lc(tess_dir, ra=ra, dec=dec, ztf=ztf)

    if ztf:
        ccd_dir = tess_dir+'s00{}/s00{}-lc-ZTF/'.format(sector, sector)
    else:
        ccd_dir = tess_dir+'s00{}/s00{}-lc/'.format(sector, sector)
        
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

    freqs, power = LombScargle(t, y).autopower(samples_per_peak=10)
    power = (power - np.mean(power)) / np.std(power)
    
    plt.figure(figsize=figsize)
    # periods = [745.87715, 704.24127, 702.89501, 701.75412,  649.39711]
    # periods = [702.89501, 351.46205]
    # labels = ['702.89501 s', '702.89501/2 s']
    periods = [701.8933, 436.9813,  352.938,  233.264]
    colors=['b', 'm', 'r', 'g', 'c']
    for i in range(len(periods)):
        # plt.axvline((1440*60.)/periods[i], c=colors[i], ls='--', lw=1, label=labels[i])
        plt.axvline((1440*60.)/periods[i], c=colors[i], ls='--', lw=1, label=str(periods[i])+' s')
        
    plt.plot(freqs, power, '-k', lw=1)
    # plt.plot(freqs, power, '.k', ms=1)

    plt.xlim([100,400])
    plt.xlabel('Frequency [1/days]')
    plt.ylabel('LS Power\nstd above mean')
    plt.legend(loc='upper center')
    plt.tight_layout()
    plt.savefig(out_dir+suffix+'_LS_Periodogram.png',dpi=300)
    print('Saved '+out_dir+suffix+'_LS_Periodogram.png')
