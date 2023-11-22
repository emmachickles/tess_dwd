import lc_utils as lcu
from Period_Finding import BLS, LS
import numpy as np
import matplotlib.pyplot as plt
import os

qmin = 0.01
# qmax = 0.15
qmax = 0.3
pmin = 2 # minutes
pmax = 10 # days
output_dir = "/home/echickle/out/"
# data_dir = "/matchfiles/data2/ATLAS/"
data_dir = '/data/ATLAS/WDs/'
wd_main = "/data/GaiaEDR3_WD_main.fits"
rp_ext = "/data/GaiaEDR3_WD_RPM_ext.fits"

gid = 6362060327430712448
fname_atlas = data_dir + str(gid)
suffix="_"+str(gid)

t, y, dy, ra, dec = lcu.load_atlas_lc(fname_atlas, pos_iqr=10, neg_iqr=10, clip=False)
print(ra)
print(dec)


plt.figure()
plt.errorbar(t, y, dy, ls=' ', c='k',capsize=1, elinewidth=1)
plt.ylabel('Flux [uJy]')
plt.xlabel('BJD')
plt.savefig(output_dir+'GID'+str(gid)+'.png', dpi=300)
print(output_dir+'GID'+str(gid)+'.png')

freqs_to_remove = []
df = 0.05
# freqs_to_remove.append([1 - df, 1 + df])
# freqs_to_remove.append([1/2. - df, 1/2. + df])
# freqs_to_remove.append([1/4. - df, 1/4. + df])

# import time
# start_time = time.time()
t, y, dy, period, bls_power_best, freqs, power, q, phi0 = \
    BLS(t,y,dy,pmin=pmin,pmax=pmax,qmin=qmin,qmax=qmax,freqs_to_remove=freqs_to_remove)
res = lcu.vet_plot(t, y, freqs, power, q, phi0, output_dir=output_dir,
                   objid=gid, objid_type='GAIAID',
                   dy=dy, suffix=suffix, wd_main=wd_main, rp_ext=rp_ext)
sig, snr, wid, period, period_min, q, phi0, epo, rp, nt, dphi = res
print(epo)
print(period*1440)
print(period*q*1440)

plt.ion()
plt.figure(figsize=(8,6))
plt.xlabel('Period [day]')
plt.ylabel('BLS Power')
plt.plot(1/freqs, power, '.k', ms=1)
plt.savefig(output_dir+'GID'+str(gid)+'_BLS_periodogram_P.png')
print(output_dir+'GID'+str(gid)+'_BLS_periodogram_P.png')
inds =  

# _, _, _, ls_period, ls_freqs, ls_power = LS(t,y,dy,freqs_to_remove=freqs_to_remove)
# suffix='_ra_{}_dec_{}'.format(ra, dec)

# res=lcu.vet_plot(t, y, ls_freqs, ls_power,output_dir=output_dir,
#                  objid=gid, objid_type=objid_type, suffix=suffix, ra=ra,
#                  dec=dec, wd_main=wd_main, rp_ext=rp_ext, 
#                  bls=False)
    
# per, q, epo = res[3], res[5], res[7]
# lcu.plot_eclipse_timing(t, y, per, epo, q, output_dir+'GAIAID_{}_{}_{}_'.format(gid, ra, dec))

from astropy.timeseries import LombScargle
# freqLS, powLS = LombScargle(t,y).autopower()
freqLS=freqs
perLS=1/freqLS
powLS=LombScargle(t,y).power(freqLS)
plt.figure(figsize=(8,6))
plt.xlabel('Period [day]')
plt.ylabel('LS Power')
plt.plot(1/freqLS, powLS, '.k', ms=1)
plt.savefig(output_dir+'GID'+str(gid)+'_LS_periodogram.png')
print(output_dir+'GID'+str(gid)+'_LS_periodogram.png')

window = np.nonzero( (perLS > 1.3) * (perLS < 1.36) )
max_ind = np.argmax(powLS[window])
print(perLS[window][max_ind])

window = np.nonzero( (perLS > 0.55) * (perLS < 0.62) )
max_ind = np.argmax(powLS[window])
print(perLS[window][max_ind])

