import numpy as np
import os
import matplotlib.pyplot as plt


mydir = '/home/echickle/data/EOF_IQR1010/'
vet_dir = mydir + 'vet_res/'
output_dir = mydir + 'met_plot/'
os.makedirs(output_dir, exist_ok=True)
suffix = '1010'

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

f = np.loadtxt(mydir+'JVR'+suffix+'.result', usecols=(1,2,3,10,11))
pow, snr, dphi = f[:,0], f[:,1], f[:,4]
wid, nt = np.int64(f[:,2]), np.int64(f[:,3])
match = np.loadtxt(mydir+'JVR'+suffix+'.result', usecols=(12), dtype='str')
inds = np.nonzero(match == 'T')
ax1.plot(pow[inds], snr[inds], '<b', label='JVR')
ax2.plot(pow[inds], wid[inds], '<b', label='JVR')
ax3.plot(nt[inds], dphi[inds], '<b', label='JVR')


if suffix == '1010':
    f = np.loadtxt(mydir+'KB'+suffix+'.result', usecols=(1,2,3,10,11))
    match = np.loadtxt(mydir+'KB'+suffix+'.result', usecols=(12), dtype='str')
if suffix == '310':
    f = np.loadtxt(mydir+'KB'+suffix+'.result', usecols=(1,2,3,11,12))
    match = np.loadtxt(mydir+'KB'+suffix+'.result', usecols=(13), dtype='str')        
pow, snr, dphi = f[:,0], f[:,1], f[:,4]
wid, nt = np.int64(f[:,2]), np.int64(f[:,3])
inds = np.nonzero(match == 'T')
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


    


