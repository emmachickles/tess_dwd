# ------------------------------------------------------------------------------

atlas_dir = "/data/ATLAS/WDs/"
output_dir='/home/echickle/out/'

# ra =102.88891
dec = 28.73981
period =   765.2038361137168 / 1440. / 60.

ra=141.06274385670
dec=35.28094464811
period = 215.2 /1440./60.
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
t, y, dy, ra, dec = lcu.load_atlas_lc(fname_atlas, clip=True, pos_iqr=5, neg_iqr=5, skiprows=0)

# ------------------------------------------------------------------------------

suffix = str(ra)+'_'+str(dec)

fig, ax = plt.subplots()
folded_t, folded_y, folded_dy = lcu.bin_timeseries(t%period, y, 75, dy=dy)
lcu.plot_phase_curve(ax, folded_t, folded_y, folded_dy, period=period,
                 ylabel="ATLAS Relative Flux")
fig.savefig(output_dir+suffix+'_binned.png', dpi=300)
print(output_dir+suffix+'_binned.png')

fig, ax=plt.subplots()
lcu.plot_phase_curve(ax, t%period, y, dy, period=period,
                     ylabel="ATLAS Relative Flux", alpha=0.6)
w = max(1, int(0.05*len(y)))
inds = np.argsort(t%period)
tconv = (t%period)[inds]
tconv = np.append(tconv, tconv+period) 
yconv = np.concatenate((y[inds], y[inds]))
tconv = np.convolve(tconv, np.ones(w), 'valid') / w
yconv = np.convolve(yconv, np.ones(w), 'valid') / w        
ax.plot(tconv*1440., yconv, '-b', lw=1)
ax.set_ylim([np.min(y)-np.std(y), np.max(y)+np.std(y)])        
fig.savefig(output_dir+suffix+'_unbinned.png', dpi=300)
print(output_dir+suffix+'_unbinned.png')

fig, ax  =plt.subplots(figsize=(10,4))
ax.errorbar(t, y, dy, fmt='.k', ms=2, elinewidth=1, capsize=1., alpha=0.8)
ax.set_xlabel('MBJD')
ax.set_ylabel('ATLAS Relative Flux')
fig.savefig(output_dir+suffix+'_raw.png', dpi=300)
print(output_dir+suffix+'_raw.png')
