import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import os
import pdb

mydir = '/home/echickle/data/EOF_TESS/'
fnames = os.listdir(mydir)
sig = []
for f in fnames:
    data = np.loadtxt(mydir+f, delimiter=',')
    if len(data.shape) == 2:
        sig.extend(data[:,3])

sig = np.array(sig)
inds = np.nonzero(~np.isnan(sig))
sig = sig[inds]
inds = np.nonzero( sig != np.inf)
sig = sig[inds]
inds = np.nonzero( sig < 150)
sig = sig[inds]

        
fig,ax = plt.subplots(figsize=(6,4))
fig.patch.set_facecolor('black')
ax.set_facecolor('black')  # Set background color
ax.spines['bottom'].set_color('w')
ax.spines['top'].set_color('w') 
ax.spines['right'].set_color('w')
ax.spines['left'].set_color('w')

plt.hist(sig, color='c', bins=50)

ax.set_xlabel("Statistical Significance of BLS Peak\n(amplitude - median) / standard deviation", color='white')  # Set xlabel color
ax.set_ylabel("White dwarf systems in TESS Sectors 56-65", color='white')  # Set ylabel color
ax.xaxis.label.set_color('white')  # Set x-axis label color
ax.yaxis.label.set_color('white')  # Set x-axis label color
ax.tick_params(axis='both', colors='white')  # Set x-axis tick label color


plt.tight_layout()


plt.savefig('/home/echickle/out/sig_hist.png', dpi=300)
print('Saved /home/echickle/out/sig_hist.png')
