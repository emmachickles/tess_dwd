# -- inputs --------------------------------------------------------------------

output_dir = '/data/submit/tess/echickle/'
img_dir    = output_dir + 's0057-png-230103-02/'
data_dir   = '/data/submit/tess/echickle/s0057-lc/'
diag_dir   = output_dir + 's0057-diag-230103-02/'

n_std = 2
wind  = 0.05

# ------------------------------------------------------------------------------

# DWD: ./sector-05/tess2018319095959-s0005-0000000471013547-0125-s_lc.fits
# https://arxiv.org/pdf/2106.15104.pdf

import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt

import pdb
import os
import sys
import gc
import fnmatch
sys.path.insert(0, "/home/submit/echickle/work/")

from KBB_Utils.KBB_Utils.Period_Finding import BLS

from astropy.io import fits
from astropy.timeseries import TimeSeries
from astropy.time import Time
from astropy.coordinates import SkyCoord
import astropy.units as u
from wotan import flatten
import lc_utils as lcu
import gc

cam = sys.argv[1]

if not os.path.exists(img_dir):
    os.makedirs(img_dir)

fail_txt = output_dir + 'cam'+str(cam)+'-failed.txt'
with open(fail_txt, 'w') as f:
    f.write('')

fnames = fnmatch.filter(os.listdir(data_dir), 'lc-'+str(cam)+'*')
    
for fname in fnames: 

    ccd = int(fname.split('-')[2][:1])
    flux = np.load(data_dir+fname)
    
    suffix = '-'+str(cam)+'-'+str(ccd)+'.npy'
    time = np.load(data_dir+'ts'+suffix)     
    coord = np.load(data_dir+'co'+suffix)
    ticid = np.load(data_dir+'id'+suffix).astype('int')        
    
    inds = np.argsort(time)
    time, flux = time[inds], flux[:,inds]

    # if cam == '2':
    #     ind = np.nonzero(ticid == 0)
    #     ticid = ticid[ind]
    #     coord = coord[ind]
    #     flux = flux[ind]
    # if cam == '4':
    #     ind = np.nonzero(ticid == 1201247611)
    #     ticid = ticid[ind]
    #     coord = coord[ind]
    #     flux = flux[ind]

    # !!
    # ind = np.nonzero(ticid == 0) # >> Kevin's DWD
    # coord = coord[ind]
    # flux = flux[ind]

    # # >> remove completed
    # fnames = os.listdir(img_dir)
    # fnames_ccd = fnmatch.filter(fnames, '*_{}_{}_*'.format(cam, ccd))
    # ticid_ccd = [int(f.split('_')[1][3:]) for f in fnames_ccd]
    # ticid_ccd = np.array(ticid_ccd)
    # inter, comm1, comm2 = np.intersect1d(ticid, ticid_ccd, return_indices=True)
    # coord = np.delete(coord, comm1, axis=0)
    # flux = np.delete(flux, comm1, axis=0)
    # ticid = np.delete(ticid, comm1) 

    # >> compute BLS
    for i in range(len(flux)):
        print(i)
        y = flux[i]
        t = time
        
        # -- prep light curve --------------------------------------------------
        t, y, flag = lcu.prep_lc(t, y, n_std=n_std, wind=0.1, diag=True,
                                 ticid=ticid[i], cam=cam, ccd=ccd,
                                 coord=coord[i], output_dir=diag_dir)
        dy = np.ones(y.shape)

        if flag:
            with open(fail_txt, 'a') as f:
                f.write(str(ticid[i])+'\n')
        else:

            # -- compute BLS ---------------------------------------------------
            t, y, dy, period, bls_power_best, freqs, power, dur, epo = \
                BLS(t,y,dy,pmin=7,pmax=0.25,qmin=0.005,qmax=0.2,remove=True)
            # if np.isnan(bls_power_best): pdb.set_trace()
            gc.collect()

            # -- plot phase curve ----------------------------------------------
            if len(fname.split('-')) > 3:
                N_ap = fname.split('-')[3]
                N_in = 'in'+fname.split('-')[4][4:7]
                N_out = 'out'+fname.split('-')[4][9:12]

                prefix = 'pow_'+str(bls_power_best)+'_per_'+str(round(period*1440,2))+\
                    '_TIC%016d'%ticid[i]+'_cam_'+str(cam)+'_ccd_'+str(ccd)+\
                    '_dur_'+str(dur)+'_epo_'+str(epo)+\
                    '_Nap_'+N_ap+'_Nin_'+N_in+'_Not_'+N_out+\
                    '_ra_{}_dec_{}_'.format(coord[i][0], coord[i][1])
            else:
                prefix = 'pow_'+str(bls_power_best)+'_per_'+str(round(period*1440,2))+\
                    '_TIC%016d'%ticid[i]+'_cam_'+str(cam)+'_ccd_'+str(ccd)+\
                    '_dur_'+str(dur)+'_epo_'+str(epo)+\
                    '_ra_{}_dec_{}_'.format(coord[i][0], coord[i][1])                

            lcu.make_phase_curve(t, y, period, dy=dy, output_dir=img_dir,
                                 prefix=prefix, freqs=freqs, power=power,
                                 ticid=ticid[i], bins=100)

                # pdb.set_trace()

# lcfile = '/scratch/submit/tess/data/tesscurl_sector_5_lc/tess2018319095959-s0005-0000000471013547-0125-s_lc.fits'

# lcfile = '/scratch/echickle/tmp/tess2021204101404-s0041-0000001201247611-0212-s_lc.fits'
# # lcfile = '/data/submit/echickle/tess2021204101404-s0041-0000001201247611-0212-s_lc.fits'
# hdul = fits.open(lcfile)
# t = hdul[1].data['TIME']
# y = hdul[1].data['PDCSAP_FLUX']

# # >> remove quality flags
# inds = np.nonzero(hdul[1].data['QUALITY'])
# t = np.delete(t, inds)
# y = np.delete(y, inds)

# inds = np.nonzero(~np.isnan(y))
# t, y = t[inds], y[inds]
# dy = np.ones(y.shape)

# # >> detrend
# # * window length: 2-3 times longer than transit duration
# window_length = 0.3
# method='biweight'
# # method='mean'
# #method=None
# if type(method) != type(None):
#     y = flatten(t, y, window_length=window_length, method=method)

# # plt.figure()
# # plt.plot(t, y, '.k', ms=1)
# # plt.xlabel('TIME')
# # plt.ylabel('PDCSAP_FLUX')
# # plt.title('TIC 1201247611, S41-20s')
# # plt.savefig('/scratch/echickle/tmp/TIC1201247611_lc.png')

# # >> pmin: minimum period (minutes)
# # >> pmax: maximum period (days)
# # >> qmin: minimum transit duraction (t_trans / period)

# pmin_array = [3]
# pmax_array = [0.1]
# qmin_array = [0.005]
# qmax_array = [0.05]
# remove_array = [True, False]

# n_iter = 0
# for pmin in pmin_array:
#     for pmax in pmax_array:
#         for qmin in qmin_array:
#             for qmax in qmax_array:
#                 for remove in remove_array:

#                     # pmin=3
#                     # # pmax=True
#                     # pmax=0.05 # = 72 minutes
#                     # # qmin = 0.005
#                     # # qmax = 0.01
#                     # qmin=0.005
#                     # qmax=0.05
#                     # remove=False
#                     t, y, dy, period, bls_power_best, freqs, power = BLS(t,y,dy,pmin=pmin,pmax=pmax,
#                                                                          qmin=qmin,qmax=qmax,remove=remove)



#                     # >> compute phase curve
#                     ts = TimeSeries(time=Time(t, format='jd'), data={'flux': y})
#                     ts_folded = ts.fold(period=period*u.d) 

#                     sorted_inds = np.argsort(ts_folded['time'])
#                     folded_t = ts_folded['time'][sorted_inds]
#                     folded_y = ts_folded['flux'][sorted_inds]

#                     fig, ax = plt.subplots(nrows=2)
#                     ax[0].set_title('TIC 1201247611, pmin: '+str(qmin)+', pmax: '+str(pmax)+'\nqmin: '+str(qmin)+', qmax: '+str(qmax)+', remove: '+str(remove)+'\nwindow_length: '+str(window_length)+', method: '+str(method)+'\nperiod: '+str(round(period*1440,2))+' min')
#                     ax[0].plot(freqs, power, '.k', ms=1)
#                     ax[0].set_xlim([0,50])
#                     ax[0].set_xlabel('Frequency')
#                     ax[0].set_ylabel('BLS Power')
#                     ax[1].plot(folded_t.value, folded_y.value, '.k', ms=1)
#                     ax[1].set_xlabel('Phase')
#                     ax[1].set_ylabel('PDCSAP_FLUX')
#                     # fname = '/data/submit/echickle/out/TIC1201247611_pc_qmin'+str(qmin)+'_qmax'+str(qmax)+'_period'+str(period)+'_sig'+str(bls_power_best)+'.png'
#                     # fname = '/data/submit/echickle/out/TIC1201247611-%02d'%n_iter+'.png'
#                     fname = img_dir+'TIC1201247611-%02d'%n_iter+'-p'+str(period)+'.png'
#                     fig.tight_layout()
#                     fig.savefig(fname)
#                     n_iter += 1

# plt.figure()
# plt.plot(freqs, power, '.k', ms=1)
# plt.xlim([0, 50])
# plt.xlabel('Frequency')
# plt.ylabel('BLS Power')
# plt.savefig('/scratch/echickle/tmp/TIC1201247611_bls.png')

# plt.figure()
# plt.plot(folded_t.value, folded_y.value, '.k', ms=1)
# plt.xlabel('Phase')
# plt.ylabel('PDCSAP_FLUX')
# plt.title('TIC 1201247611, pmin: '+str(qmin)+', pmax: '+str(pmax)+', qmin: '+str(qmin)+', qmax: '+str(qmax))
# fname = '/scratch/echickle/tmp/TIC1201247611_pc_qmin'+str(qmin)+'_qmax'+str(qmax)+'_period'+str(period)+'_sig'+str(bls_power_best)+'.png'
# plt.savefig(fname)
# print('Saved '+fname)

# with open('/home/submit/echickle/foo.txt', 'w') as f:
#     f.write('Period '+str(period)+'\nPower '+str(bls_power_best))
# >> load white dwarf catalog
# wd_cat  = pd.read_csv('/home/echickle/work/tess_dwd/WDs.txt', header=None, sep='\s+')

# >> WDs observed in S30 short cadence
# s30 = np.loadtxt('/home/echickle/work/tess_dwd/all_targets_S030_v1.txt')
# n = len(np.intersect1d(s30[:,0], wd_cat[0]))



# >> sigma clipping exercise
# for i in [5, 7, 8, 10]:
#     npts = []
#     for j in range(len(flux)):
#         cond1 = np.count_nonzero(flux[j] < np.mean(flux[j]) - i*np.std(flux[j]))
#         cond2 = np.count_nonzero(flux[j] > np.mean(flux[j]) + i*np.std(flux[j]))        
#         npts.append( cond1+cond2 )
#     plt.figure()
#     plt.hist(npts, bins=100)
#     plt.xlabel('Number of points clipped (sigma='+str(i)+')')
#     plt.ylabel('Number of targets')
#     plt.savefig('/data/submit/echickle/out/imgs4/sigma_'+str(i)+'.png')
# pdb.set_trace()

# inds = np.nonzero(ticid == 1201247611)
# flux = flux[inds]
# ticid = ticid[inds]


# for i in range(len(flux)):
#     plt.figure()
#     plt.plot(time, flux[i], '.k')
#     plt.title('TIC '+str(ticid[i]))
#     plt.xlabel('Time')
#     plt.ylabel('Flux')
#     plt.savefig(output_dir+'sector_41_raw_lcs/lc_'+str(ticid[i])+'.png')
#     plt.close()
#     print(i)
# # <<


# pdb.set_trace()


# fnames = []
# for i in os.listdir(data_dir):
#     for j in os.listdir(data_dir+i):
#         fnames.append(data_dir+i+'/'+j)

# with open(output_dir+'bls_power.txt', 'w') as f:
#     f.write('Sector,TICID,period,bls_power_best\n')

# # for fname in [data_dir + 'sector-05/tess2018319095959-s0005-0000000471013547-0125-s_lc.fits']:
# for fname in fnames:
#     try:
#         with fits.open(fname, memmap=False) as hdul:
#             y = hdul[1].data['PDCSAP_FLUX'] # >> unit = e-/s
#             dy = hdul[1].data['PDCSAP_FLUX_ERR'] # >> unit = e-/s
#             t = hdul[1].data['TIME'] # >> unit = BJD - 2457000, days

#     except:
#         print('Failed to open the following FITS file:')
#         print(fname)
#     gc.collect()

#     # >> mask out nans
#     inds = np.nonzero(~np.isnan(y))
#     y, dy, t = y[inds], dy[inds], t[inds]

#     # >> compute BLS
#     t, y, dy, period, bls_power_best = BLS(t,y,dy,pmin=3,pmax=0.1,qmin=0.05,qmax=0.1,remove=True)
#     # >> pmin: minimum period in minutes, default 3 minutes
#     # >> pmax: maximum period in days, default 4/baseline

#     s = int(fname.split('-')[2][1:])
#     tic = int(fname.split('-')[3])
#     with open(output_dir+'bls_power.txt', 'a') as f:
#         f.write(str(s)+','+str(tic)+','+str(period)+','+str(bls_power_best)+'\n')


#     # >> compute phase curve
#     ts = TimeSeries(time=Time(t, format='jd'), data={'flux': y})
#     ts_folded = ts.fold(period=period*u.d) 

#     sorted_inds = np.argsort(ts_folded['time'])
#     folded_t = ts_folded['time'][sorted_inds]
#     folded_y = ts_folded['flux'][sorted_inds]

#     plt.figure()
#     plt.plot(folded_t.value, folded_y.value, '.k', ms=1)
#     plt.xlabel('Phase')
#     plt.ylabel('Flux')
#     plt.savefig(img_dir+str(bls_power_best)+'-'+fname.split('-')[2]+'-'+fname.split('-')[3]+'-'+str(period)+'_folded.png')
            # -- check if failure --------------------------------------------------
            # if not os.path.exists(img_dir+prefix+'phase_curve.png'):
            #     with open(fail_txt, 'a') as f:
            #         f.write('TIC'+str(ticid[i])+'\n')
        # # >> remove nans
        # inds = np.nonzero(~np.isnan(y))
        # t, y = t[inds], y[inds]

        # if y.shape[0] < 1000:
        # else:


        #     med = np.median(y)
        #     std = np.std(y)

        #     # >> sigma-clip to 5 sigma
        #     inds = np.nonzero( (y > med - 5*std) * (y < med + 5*std) )
        #     t, y = t[inds], y[inds]

        #     # >> normalize        
        #     if np.median(y) < 0:
        #         y = y / np.abs(med) + 2.
        #     else:
        #         y = y/med

        #     # >> detrending 
        #     y = flatten(t, y, window_length=0.1, method='biweight') # >> can't handle negative flux
        #     # >> remove
        #     if np.count_nonzero(np.isnan(y)) < 2000:
        #         inds = np.nonzero(~np.isnan(y))
        #         t, y = t[inds], y[inds]


