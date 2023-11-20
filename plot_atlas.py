# ------------------------------------------------------------------------------

atlas_dir = "/data/ATLAS/WDs/"
output_dir='/home/echickle/out/'

# ra =102.88891
dec = 28.73981
period =   765.2038361137168 / 1440. / 60.

ra=141.06274385670
dec=35.28094464811
period = 215.2 /1440./60.

ra = 	240.8656161352
dec = -65.8640973300
# period = 0.13206648882803174
period = 190.17308525 / 1440

ra = 215.6801
dec = -44.07953
period = 67.95327296/1440

ra =273.2963257
dec =42.86401367
period=0.03552830646

# ra = 44.3044
# dec = -40.65052
# period = 41.68733132/1440

ra = 44.3044
dec = -40.65052
period = 41.68733132/1440

ra = 221.25516830820
dec = -36.56763617409
period = 28.63859828*2/1440

ra = 311.14316
dec = -78.70051
period = 1.336833166965929
# period = 0.571572277476107
# period = 0.3506

bins=100
# figsize=(5,2)
figsize=(4,3)

# ------------------------------------------------------------------------------

import lc_utils as lcu
import numpy as np
import matplotlib.pyplot as plt
import sys
import os
import pdb
os.makedirs(output_dir, exist_ok=True)

# ------------------------------------------------------------------------------

# load data
fname_atlas = lcu.get_atlas_lc(atlas_dir, ra=ra, dec=dec)
t, y, dy, ra, dec = lcu.load_atlas_lc(fname_atlas, clip=True, pos_iqr=3, neg_iqr=10, skiprows=0)

# ------------------------------------------------------------------------------

suffix = str(ra)+'_'+str(dec)

import matplotlib as mpl
mpl.rcParams['font.family'] = 'serif'

fig, ax = plt.subplots(figsize=figsize)

t0 = 38/1440
t = t - t0

folded_t, folded_y, folded_dy = lcu.bin_timeseries(t%period, y, bins, dy=dy)
#  folded_t = folded_t # - 7/1440 # !!

np.savetxt(output_dir+suffix+'_data.txt', np.array([t,y,dy]).T)

folded_y = folded_y - 2.5*np.min(folded_y) # !!
folded_dy = folded_dy  /np.median(folded_y)
folded_y = folded_y / np.median(folded_y)
# lcu.plot_phase_curve(ax, folded_t, folded_y, folded_dy, period=period,
#                  ylabel="ATLAS Relative Flux")
lcu.plot_phase_curve(ax, folded_t, folded_y, folded_dy, err_sigma=5,
                 ylabel="Flux (normalised)")
#                  ylabel="ATLAS Relative Flux")

# ax.errorbar(folded_t, folded_y, folded_dy, ls=' ', capsize=1, fmt='.k',
#             elinewidth=1, ms=2)
# ax.set_ylabel('ATLAS Relative Flux')
plt.tight_layout()
# fig.savefig(output_dir+suffix+'_binned.png', dpi=300)
# print(output_dir+suffix+'_binned.png')
fig.savefig(output_dir+suffix+'_binned_P='+str(period)+'.png', dpi=300)
print(output_dir+suffix+'_binned_P='+str(period)+'.png')


fig, ax=plt.subplots(figsize=figsize)
lcu.plot_phase_curve(ax, t%period, y, dy, period=period,
                     ylabel="ATLAS Relative Flux", alpha=0.6,
                     err_sigma=5)
w = max(1, int(0.05*len(y)))
inds = np.argsort(t%period)
tconv = (t%period)[inds]
tconv = np.append(tconv, tconv+period) 
yconv = np.concatenate((y[inds], y[inds]))
tconv = np.convolve(tconv, np.ones(w), 'valid') / w
yconv = np.convolve(yconv, np.ones(w), 'valid') / w        
ax.plot(tconv*1440., yconv, '-b', lw=1)
ax.set_ylim([np.min(y)-np.std(y), np.max(y)+np.std(y)])        
# fig.savefig(output_dir+suffix+'_unbinned.png', dpi=300)
# print(output_dir+suffix+'_unbinned.png')

fig.savefig(output_dir+suffix+'_unbinned_P='+str(period)+'.png', dpi=300)
print(output_dir+suffix+'_unbinned_P='+str(period)+'.png')


fig, ax  =plt.subplots(figsize=figsize)
ax.errorbar(t, y, dy, fmt='.k', ms=2, elinewidth=1, capsize=1., alpha=0.8)
ax.set_xlabel('MBJD')
ax.set_ylabel('ATLAS Relative Flux')
fig.savefig(output_dir+suffix+'_raw.png', dpi=300)
print(output_dir+suffix+'_raw.png')

