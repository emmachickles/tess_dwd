import os
import pdb
import numpy as np
import matplotlib.pyplot as plt
from wotan import flatten

data_dir = '/home/echickle/data/s0056/s0056-lc-bkgsub/'
qflag_dir = '/home/echickle/data/QLPqflags/sec56/'
suffix = '-2-2.npy'
out_dir = '/home/echickle/out/qflag/'
os.makedirs(out_dir, exist_ok=True)
out1 = '/home/echickle/out/processing/'
os.makedirs(out1, exist_ok=True)

n_std = 5
detrend = "wotan"
wind = 0.1

t = np.load(data_dir+'ts'+suffix)
cn = np.load(data_dir+'cn'+suffix)
co = np.load(data_dir+'co'+suffix)
zid = np.load(data_dir+'id'+suffix)
file_names = ['orbit119cam2ccd2_qflag.txt', 'orbit120cam2ccd2_qflag.txt']
qflag_data = []
for f in file_names:
    qflag_data.extend(np.loadtxt(qflag_dir+f))
qflag_data = np.array(qflag_data)
bad_inds = np.nonzero(qflag_data[:,1])[0]
bad_cadence = qflag_data[:,0][bad_inds]
_, comm1, comm2 = np.intersect1d(cn, bad_cadence, return_indices=True)


file_names = [f[7:-4] for f in os.listdir(data_dir) if 'lc-2-2-' in f]
# file_names = [0]
for f in file_names:
    y_list = np.load(data_dir+'lc-2-2-'+f+'.npy')
    # y_list = np.load(data_dir+'lc'+suffix)
    for i in range(len(y_list)):
        y = y_list[i]

        plt.figure()
        plt.plot(t, y, '.k', ms=1)
        plt.plot(t[comm1], y[comm1], 'xr', ms=3)
        plt.title(zid[i])
        plt.xlabel('Time [BTJD]')
        plt.ylabel('Background-subtracted flux [electrons/s]')
        plt.savefig(out_dir+'phot_bkgsub_{}_{}_{}.png'.format(co[i][0], co[i][1], f), dpi=300)
        plt.close()
        print(out_dir+'phot_bkgsub_{}_{}_{}.png'.format(co[i][0], co[i][1], f))
        
        fig, ax = plt.subplots(nrows=3, sharex=True, figsize=(8,8))
        ax[-1].set_xlabel('Time [BTJD]')
        
        t1, y1 = np.delete(t, comm1), np.delete(y, comm1) # Remove non-zero quality flags
        
        inds = np.argsort(t1)
        t1, y1 = t1[inds], y1[inds] # Sort time array

        med = np.median(y1)        
        y1 = y1 / med # Normalize light curve
        if med < 0:
            y1 = -1 * y1 + 2
        if np.min(y1) < 0:
            y1 -= np.min(y1)
        ax[0].plot(t1, y1, '.k', ms=1)
        ax[0].set_ylabel('Relative flux')
        
        y1 = flatten(t1, y1, window_length=wind, method='biweight')
        inds = np.nonzero(~np.isnan(y1))
        t1, y1 = t1[inds], y1[inds]
        if len(t1) == 0:
            pdb.set_trace()
        ax[1].plot(t1, y1, '.k', ms=1)
        ax[1].set_ylabel('Detrended flux')        

        med = np.median(y1)
        std = np.std(y1)
        inds = np.nonzero( (y1 > med - n_std*std) * (y1 < med + n_std*std) )
        t1, y1 = t1[inds], y1[inds]
        ax[2].plot(t1, y1, '.k', ms=1)
        ax[2].set_ylabel('Sigma-clipped flux')

        plt.tight_layout()
        fig.savefig(out1+'phot_bkgsub_{}_{}_{}.png'.format(co[i][0], co[i][1], f), dpi=300)
        plt.close()
        print(out1+'phot_bkgsub_{}_{}_{}.png'.format(co[i][0], co[i][1], f))
        
