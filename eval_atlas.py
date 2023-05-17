import numpy as np
import os
import matplotlib.pyplot as plt


mydir = '/home/echickle/data/EOF_IQR1010/'
vet_dir = mydir + 'vet_res/'
output_dir = mydir + 'met_plot/'
os.makedirs(output_dir, exist_ok=True)

pow_all, snr_all, wid_all, nt_all, dphi_all = [], [], [], [], []
fnames = os.listdir(vet_dir)
for i in range(len(fnames)):
    # gid, pow, snr, wid, per, q, phi0, nt, dphi
    f = np.loadtxt(vet_dir + fnames[i])
    pow_all.extend(f[:,1])
    snr_all.extend(f[:,2])
    wid_all.extend(np.int64(f[:,3]))
    nt_all.extend(np.int64(f[:,7]))
    dphi_all.extend(f[:,8])

fig1, ax1 = plt.subplots()
ax1.plot(pow_all, snr_all, '.k', ms=1, alpha=0.5)
ax1.set_xlabel('Peak Significance = (peak - median)/MAD')
ax1.set_ylabel('SNR = depth/MAD')
fig2, ax2 = plt.subplots()
ax2.plot(pow_all, wid_all, '.k', ms=1, alpha=0.5)
ax2.set_xlabel('Peak Significance = (peak - median)/MAD')
ax2.set_ylabel('Peak width')
fig3, ax3 = plt.subplots()
ax3.plot(nt_all, dphi_all, '.k', ms=1, alpha=0.5)
ax3.set_xlabel('Number of points in-eclipse')
ax3.set_ylabel('Max(diff(phi))')

# f = np.loadtxt(mydir+'JVR'+suffix+'.result', usecols=(1,2,3,10,11))
# pow, snr, dphi = f[:,0], f[:,1], f[:,4]
# wid, nt = np.int64(f[:,2]), np.int64(f[:,3])
# match = np.loadtxt(mydir+'JVR'+suffix+'.result', usecols=(12), dtype='str')
# inds = np.nonzero(match == 'T')
# ax1.plot(pow[inds], snr[inds], '<b', label='JVR')
# ax2.plot(pow[inds], wid[inds], '<b', label='JVR')
# ax3.plot(nt[inds], dphi[inds], '<b', label='JVR')



# f cols: ra dec sig snr wid per q phi0 nt dphi
f = np.loadtxt('GPU_UCB.result', skiprows=1, usecols=(1,2,3,4,5,6,8,9,12,13))
pow, snr, dphi = f[:,2], f[:,3], f[:,9]
wid, nt = np.int64(f[:,4]), np.int64(f[:,8])
per = f[:,5]

cat = np.loadtxt("Kevin's UCBs - UCBs.csv", delimiter=',', skiprows=1, usecols=(2,3,4))
from astropy.coordinates import SkyCoord
import astropy.units as u
co = SkyCoord(ra=f[:,0], dec=f[:,1], unit=(u.degree, u.degree), frame='icrs')
cat = SkyCoord(ra=cat[:,0], dec=cat[:,1], unit=(u.degree, u.degree), frame='icrs')
idx, d2d, _ = co.match_to_catalog_sky(cat)
good_idx = idx[np.nonzero(d2d.to(u.arcsec).value < 2)]

match = []
no_match = []
per_true = cat[:,2][good_idx]
import pdb
pdb.set_trace()
dt = 2/1440
for i in range(len(f)):
    if (per > per_true[i]-dt and per < per_true[i]+dt) or\
    (per > per_true[i]*2-dt and per < per_true[i]*2+dt) or\
    (per > per_true[i]/2-dt and per < per_true[i]/2+dt):
        match.append(i)
    else:
        no_match.append(i)
    

ax1.plot(pow[inds], snr[inds], '>g', label='KBUCB')
ax2.plot(pow[inds], wid[inds], '>g', label='KBUCB')
ax3.plot(nt[inds], dphi[inds], '>g', label='KBUCB')

f = np.loadtxt(mydir+'KBWD'+suffix+'.result', usecols=(1,2,3,11,12), delimiter=',')
pow, snr, dphi = f[:,0], f[:,1], f[:,4]
wid, nt = np.int64(f[:,2]), np.int64(f[:,3])
match = np.loadtxt(mydir+'KBWD'+suffix+'.result', usecols=(13), delimiter=',', dtype='str')
inds = np.nonzero(match == 'T')
ax1.plot(pow[inds], snr[inds], '>g', label='KBUCB')
ax2.plot(pow[inds], wid[inds], '>g', label='KBUCB')
ax3.plot(nt[inds], dphi[inds], '>g', label='KBUCB')

f = np.loadtxt(mydir+'TESS'+suffix+'.result', usecols=(1,2,3,11,12), delimiter=',')
pow, snr, dphi = f[:,0], f[:,1], f[:,4]
wid, nt = np.int64(f[:,2]), np.int64(f[:,3])
match = np.loadtxt(mydir+'TESS'+suffix+'.result', usecols=(13), delimiter=',', dtype='str')
inds = np.nonzero(match == 'T')
ax1.plot(pow[inds], snr[inds], 'vr', label='TESS')
ax2.plot(pow[inds], wid[inds], 'vr', label='TESS')
ax3.plot(nt[inds], dphi[inds], 'vr', label='TESS')

ax1.legend()
fig1.tight_layout()
fig1.savefig(output_dir+'pow_snr_'+suffix+'.png')
print(output_dir+'pow_snr_'+suffix+'.png')
ax1.set_xlim([0,425])
ax1.set_ylim([0,15])
fig1.savefig(output_dir+'pow_snr_'+suffix+'_zoom.png')

ax2.legend()
fig2.tight_layout()
fig2.savefig(output_dir+'pow_wid_'+suffix+'.png')
print(output_dir+'pow_wid_'+suffix+'.png')
ax2.set_xlim([0,425])
ax2.set_ylim([0,175])
fig2.savefig(output_dir+'pow_wid_'+suffix+'_zoom.png')

ax3.legend()
fig3.tight_layout()
fig3.savefig(output_dir+'nt_dphi_'+suffix+'.png')
print(output_dir+'nt_dphi_'+suffix+'.png')   
ax3.set_xlim([0,600])
ax3.set_ylim([0,0.05])
fig3.savefig(output_dir+'nt_dphi_'+suffix+'_zoom.png')


    


