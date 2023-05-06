import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

sector = 57
pmax = 0.15*1440
out_dir = "/scratch/echickle/s%04d-ZTF/"%sector
cat = np.loadtxt("/scratch/echickle/ZTF_Eclipses.txt", dtype='str')
ra_ztf, dec_ztf = cat[:,1].astype('float'), cat[:,2].astype('float')
per_ztf = cat[:,3].astype('float')*1440

power = [] 
snr = []
wid = []
per = []
nt = []
dphi = []
per_ztf_list = []
ra_list = []
dec_list = []
loc = []
match = []


for cam in [1,2,3,4]:
    for ccd in [1,2,3,4]:
        ccd_dir = out_dir + 's%04d'%sector + '-bls-{}-{}/'.format(cam, ccd)

        f = '/scratch/echickle/s0057-ZTF/s0057-gpu-res/GPU-57-{}-{}.result'.format(cam, ccd)
        cat = pd.read_csv(f, header=None)
        cat_ticid = cat[0].to_numpy().astype('int')

        fnames = os.listdir(ccd_dir)
        for i in range(len(fnames)):
            fsplit = fnames[i].split('_')
            ra = float(fsplit[19])
            dec = float(fsplit[21])
            ind = np.nonzero(ra_ztf == ra)
            if per_ztf[ind] < pmax:
                p = float(fsplit[7])
                if np.abs(per_ztf[ind] - p) < 3:
                    match.append(True)
                else:
                    match.append(False)

                power.append(float(fsplit[1]))
                snr.append(float(fsplit[3]))
                wid.append(float(fsplit[5]))
                per.append(p)
                per_ztf_list.append(per_ztf[ind][0])
                ra_list.append(ra)
                dec_list.append(dec)
                loc.append('{} {} {}'.format(sector, cam, ccd))

                ticid = int(fsplit[12][3:])
                ind = np.nonzero(cat_ticid == ticid)[0][0]
                nt.append( cat.iloc[ind][11] )
                dphi.append(cat.iloc[ind][12])

                

power, snr, wid, per = np.array(power), np.array(snr), np.array(wid), np.array(per)
nt, dphi = np.array(nt), np.array(dphi)
per_ztf_list, loc = np.array(per_ztf_list), np.array(loc)
ra_list, dec_list, match = np.array(ra_list), np.array(dec_list), np.array(match)
res = np.array([ra_list, dec_list, loc, power, snr, wid, per, per_ztf_list, match], dtype='str').T                
np.savetxt(out_dir+'eval_jvr.txt', res, fmt='%s',
           header='RA Dec Pow SNR Wid Per Per_ZTF Match')

# # >> signals discovered in S61

ticid_tess = np.array([803489769, 36085812, 800042858, 270679102, 455206965, 452954413, 767706310, 96852226, 677736827])
pow_tess = []
snr_tess = []
wid_tess = []
per_tess = []
nt_tess = []
dphi_tess= []

# >> all signals (real and fake) in S61
power_all = []
snr_all = []
wid_all = []
per_all = []
nt_all = []
dphi_all = []
for cam in [1,2,3,4]:
    for ccd in [1,2,3,4]:
        f = '/scratch/echickle/s0061/s0061-gpu-res/GPU-61-{}-{}.result'.format(cam, ccd)
        cat = pd.read_csv(f, header=None)
        power_all.extend(cat[1])
        snr_all.extend(cat[2])
        wid_all.extend(cat[3])
        per_all.extend(cat[4])
        nt_all.extend(cat[11])
        dphi_all.extend(cat[12])

        cat_ticid = cat[0].to_numpy().astype('int')
        _, inds, _ = np.intersect1d(cat_ticid, ticid_tess, return_indices=True)
        pow_tess.extend(cat[1].to_numpy()[inds])
        snr_tess.extend(cat[2].to_numpy()[inds])
        wid_tess.extend(cat[3].to_numpy()[inds])
        per_tess.extend(cat[4].to_numpy()[inds])
        nt_tess.extend(cat[11].to_numpy()[inds])
        dphi_tess.extend(cat[12].to_numpy()[inds])

power_all = np.array(power_all)

plt.figure()
plt.plot(power_all, snr_all, '.k', ms=1, alpha=0.5)
plt.plot(pow_tess, snr_tess, 'vr', label='TESS')

inds = np.nonzero(match)
plt.plot(power[inds], snr[inds], '<b', label='ZTF-JVR')
plt.xlabel('Peak Significance = (peak - median)/MAD')
plt.ylabel('SNR = depth/MAD')
plt.legend()
plt.tight_layout()
plt.savefig(out_dir + 'pow_snr.png')
print(out_dir + 'pow_snr.png')

plt.xlim([0, 200])
plt.ylim([0, max(snr_tess)])
plt.savefig(out_dir + 'pow_snr_zoom.png')
print(out_dir + 'pow_snr_zoom.png')

plt.figure()
plt.plot(power_all, wid_all, '.k', ms=1, alpha=0.5)
plt.plot(pow_tess, wid_tess, 'vr', label='TESS')
plt.plot(power[inds], wid[inds], '<b', label='ZTF-JVR')
plt.xlabel('Peak Significance = (peak - median)/MAD')
plt.ylabel('Peak width')
plt.legend()
plt.tight_layout()
plt.savefig(out_dir + 'pow_wid.png')
print(out_dir + 'pow_wid.png')

plt.xlim([0, 200])
plt.ylim([0, max(wid_tess)])
plt.savefig(out_dir + 'pow_wid_zoom.png')
print(out_dir + 'pow_wid_zoom.png')


plt.figure()
plt.plot(nt_all, dphi_all, '.k', ms=1, alpha=0.5)
plt.plot(nt_tess, dphi_tess, 'vr', label='TESS')
plt.plot(nt[inds], dphi[inds], '<b', label='ZTF-JVR')
plt.xlabel('Number of points in transit')
plt.ylabel('Dphi')
plt.legend()
plt.tight_layout()
plt.savefig(out_dir + 'nt_dphi.png')
print(out_dir + 'nt_dphi.png')

plt.figure()
_=plt.hist(power_all[np.nonzero(power_all < 3000)], bins=100)
plt.xlabel('Peak Significance')
plt.ylabel('Num objects in TESS S61')
plt.savefig(out_dir + 'pow_hist.png')
print(out_dir + 'pow_hist.png')
