from astropy.io import fits
import lc_utils as lcu
import pandas as pd 
import pdb
import os
import numpy as np
import matplotlib.pyplot as plt
from wotan import flatten
import sys
sys.path.insert(0, "/home/submit/echickle/work/")    
from Period_Finding import BLS, LS_Astropy

# test

n_std=5
wind=0.05

data_dir = '/data/submit/tess/echickle/KBUCB-lc/'
out_dir = '/data/submit/tess/echickle/KBUCB-png-230308/'
os.makedirs(out_dir, exist_ok=True)

ucb_catalog = pd.read_csv('/data/submit/echickle/KB_UCBs.tsv', sep='\t', header=1, index_col=False)

fnames = os.listdir(data_dir)
ztfid_list = []
suffix_list = []
for f in fnames:
    if 'co-' in f:
        ztfid_list.append(f.split('-')[3][:-4])
        suffix_list.append(f[2:])

for i in range(len(ztfid_list)):
    ztfid, suffix = ztfid_list[i], suffix_list[i]
    print(ztfid)
    t = np.load(data_dir+'ts'+suffix)
    y = np.load(data_dir+'lc'+suffix)[0]
    coord = np.load(data_dir+'co'+suffix)[0]

    per = float(ucb_catalog.loc[ucb_catalog['ZTF Name'] == ztfid]['Period (days)'])

    t, y, flag = lcu.prep_lc(t, y, n_std=n_std, wind=wind)

    # >> BLS
    dy = np.ones(y.shape)
    # _,_,_, period, bls_power_best, freqs, power, dur, epo = \
    #     BLS(t,y,dy,pmin=7,pmax=0.25,qmin=0.1,qmax=0.2,remove=True)

    _, _, _, period, bls_power_best, freqs, power = LS_Astropy(t,y,dy,pmax=0.25)
    
    # >> fold
    t = t % per * 1440
    inds = np.argsort(t)
    t, y = t[inds], y[inds]

    # >> bin
    dy = np.ones(y.shape)
    t, y, dy = lcu.bin_timeseries(t, y, 60, dy)

    shift = np.max(t) - np.min(t)
    plt.figure(figsize=(5,3))
    plt.errorbar(t, y, yerr=np.ones(len(y))*0.1,
                   fmt='.k', ms=1, elinewidth=1)
    plt.errorbar(t+shift, y, yerr=np.ones(len(y))*0.1,
                   fmt='.k', ms=1, elinewidth=1)
    plt.xlabel('Time [minutes]')
    plt.ylabel('Relative Flux')
    plt.tight_layout()
    # plt.figure()
    # plt.plot(t, y, '.k', ms=1)
    # shift = np.max(t) - np.min(t)        
    # plt.plot(t+shift, y, '.k', ms=1) 
    # plt.xlabel('Time [minutes]')
    # plt.ylabel('Relative Flux')
    plt.savefig(out_dir + suffix[1:-4]+'_phase_curve.png', dpi=300)

    plt.figure()
    plt.plot(freqs, power, '-k', lw=0.5, alpha=0.7)
    plt.xlabel('Frequency [1/days]')
    plt.ylabel('BLS Power')
    plt.title('period: '+str(round(period*1440,2))+' min')
    plt.savefig(out_dir+suffix[1:-4]+'_spectrum.png')


    # plt.figure()
    # plt.axvline(per*1440,color='r', linestyle='dashed')        
    # plt.plot(1440/freqs, power, '-k', lw=0.5, alpha=0.7)
    # plt.xlim([400/60,120])
    # plt.xlabel('Period [minutes]')
    # plt.ylabel('BLS Power')
    # plt.savefig(out_dir+suffix[1:-4]+'_spectrum_zoom.png')


    plt.figure(figsize=(5,3))
    plt.title('period :{}'.format(np.round(period*1440,2)))
    plt.axvline(per*1440,color='r', linestyle='dashed')            
    plt.plot(1440/freqs, power, '-k', lw=0.5, alpha=0.7)
    plt.xlim([400/60,120])    
    plt.xlabel('Period [minutes]')
    plt.ylabel('BLS Power')
    plt.tight_layout()
    # plt.savefig('/home/submit/echickle/gaia14aae_tglc_spectrum_zoom.png', dpi=300)
    plt.savefig(out_dir+suffix[1:-4]+'_spectrum_zoom.png', dpi=300)    
    
    plt.figure()
    plt.axvline(per*1440,color='r', linestyle='dashed')    
    plt.plot(1440/freqs, power, '-k', lw=0.5, alpha=0.7)
    # plt.xlim([0,120])
    plt.xlim([per*1440-1, per*1440+1])    
    plt.xlabel('Period [minutes]')
    plt.ylabel('BLS Power')
    plt.savefig(out_dir+suffix[1:-4]+'_spectrum_zzoom.png')
    print(out_dir+suffix[1:-4]+'_spectrum_zzoom.png')

    plt.figure()
    plt.hist(power,bins=150)
    plt.xlabel('Power')
    plt.ylabel('Frequency bins')
    plt.xscale('log')
    plt.savefig(out_dir+suffix[1:-4]+'_hist.png')
