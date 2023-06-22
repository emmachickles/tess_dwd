import os
import pandas as pd
import numpy as np
import lc_utils as lcu
import matplotlib.pyplot as plt

sector = 56
n_std = 5
detrend = "wotan"
wind = 0.1

data_dir = '/scratch/data/tess/lcur/ffi/s%04d-lc/'%sector
output_dir = '/scratch/echickle/s%04d/'%sector
gpu_dir = output_dir + 's%04d'%sector + '-gpu-res/'
vet_dir = output_dir + 's%04d'%sector + '-vet-res/'
os.makedirs(vet_dir, exist_ok=True)

result = []

for cam in [1,2,3,4]:
    for ccd in [1,2,3,4]:
        print('Cam {} CCD {}'.format(cam, ccd))

        suffix = '-'+str(cam)+'-'+str(ccd)+'.npy'
        ticid_ccd = np.load(data_dir+'id'+suffix).astype('int')
        coord = np.load(data_dir+'co'+suffix)

        flux = np.load(data_dir+'lc'+suffix)

        gpu_f = gpu_dir + 'GPU-{}-{}-{}.result'.format(sector, cam, ccd)

        # ticid, ra, dec, sig, snr, wid, period, period_min, q, phi0, dur, epo, rp, nt, dphi
        # ticid, sig, snr, wid, period, period_min, q, phi0
        cat = pd.read_csv(gpu_f, header=None) # !!
        for i in range(len(cat)): 
            t = np.load(data_dir+'ts'+suffix)     

            # !! 
            ticid = int(cat.iloc[i][0])
            pow = cat.iloc[i][1]
            period = cat.iloc[i][4]
            q = cat.iloc[i][6]
            phi0 = cat.iloc[i][7]

            ind = np.nonzero(ticid_ccd == ticid)[0][0]
            co, y = coord[ind], flux[ind]
            t, y, flag = lcu.prep_lc(t, y, n_std=n_std, detrend=detrend, wind=wind)

            if not flag:
                snr, phi, _, _, _, nt, dphi = lcu.calc_snr(t, y, period, q, phi0)
                dphi_med = np.median(np.diff(np.sort(phi)))
                res = [ticid, co[0], co[1], pow, period, snr, nt, dphi, dphi_med]
                result.append(res)


result = np.array(result)
np.savetxt(vet_dir+'VET.result', result,
           fmt=('%10.5f,'*len(res))[:-1],
           header='ticid, ra, dec, sig, period, snr, nt, dphi, dphi_med')
print(vet_dir+'VET.result')

ticid_tess = np.array([803489769, 36085812, 800042858, 270679102, 455206965, 452954413, 767706310, 96852226, 677736827])
_, inds, _ = np.intersect1d(result[:,0].astype('int'), ticid_tess, return_indices=True)
result_tess = result[inds]

plt.figure()
plt.plot(result[:,3], result[:,5], '.k', ms=1, alpha=0.5)
plt.plot(result_tess[:,3], result_tess[:,5], 'vr', label='TESS')
plt.legend()
plt.xlabel('Power')
plt.ylabel('SNR')
plt.xlim([0, 200])
plt.savefig(vet_dir+'pow_snr.png')
print(vet_dir+'pow_snr.png')

plt.figure()
plt.plot(result[:,3], result[:,6], '.k', ms=1, alpha=0.5)
plt.plot(result_tess[:,3], result_tess[:,6], 'vr', label='TESS')
plt.legend()
plt.xlabel('Power')
plt.ylabel('nt')
plt.xlim([0, 200])
plt.savefig(vet_dir+'pow_nt.png')
print(vet_dir+'pow_nt.png')

plt.figure()
plt.plot(result[:,3], result[:,7], '.k', ms=1, alpha=0.5)
plt.plot(result_tess[:,3], result_tess[:,7], 'vr', label='TESS')
plt.legend()
plt.xlabel('Power')
plt.ylabel('Max(dphi)')
plt.xlim([0, 200])
plt.savefig(vet_dir+'pow_dphi_max.png')
print(vet_dir+'pow_dphi_max.png')

plt.figure()
plt.plot(result[:,3], result[:,8], '.k', ms=1, alpha=0.5)
plt.plot(result_tess[:,3], result_tess[:,8], 'vr', label='TESS')
plt.legend()
plt.xlabel('Power')
plt.ylabel('Median(dphi)')
plt.xlim([0, 200])
plt.tight_layout()
plt.savefig(vet_dir+'pow_dphi_med.png')
print(vet_dir+'pow_dphi_med.png')


