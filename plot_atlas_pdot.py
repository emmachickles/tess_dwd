import lc_utils as lcu
import numpy as np
import sys
sys.path.insert(0, '/home/echickle/work/KBB_Utils/KBB_Utils/')
import LC_Tools as lct
import matplotlib.pyplot as plt

atlas_dir = "/data/ATLAS/WDs/"
output_dir = '/home/echickle/'
clip = True
skiprows = 1
objid_type=None

# P = 0.004800828
# Pdot = -2.373e-11 # s s^-1
# t0 = 2458305.6827886 - 2400000.5
# fname_atlas = '/home/echickle/data/atlasforcedphotometryresults_UCB/job465518.txt'

P = 765.206543 / 1440. / 60.
Pdot = -9.8e-12
t0 = 0
# fname_atlas = '/home/echickle/data/atlasforcedphotometryresults_UCB/job465520.txt'
fname_atlas = atlas_dir + '887586156003183744'

# -- load light curve ----------------------------------------------------------

t, y, dy, ra, dec = lcu.load_atlas_lc(fname_atlas, clip=True, pos_iqr=2, neg_iqr=10,  skiprows=skiprows)
t -= 15./86400


binned_LC = lct.binning(t,y,dy,P,Pdot=Pdot,t0=t0, N=100) # phases, y, dy
# binned_LC = lct.binning(t,y,dy,P,t0=t0, N=100) # phases, y, dy

# -- plot phase curves ---------------------------------------------------------

fig, ax = plt.subplots()
ax.errorbar(binned_LC[:,0], binned_LC[:,1], binned_LC[:,2], fmt='.k', ms=2,
             elinewidth=1, capsize=1., alpha=0.5)
ax.set_xlabel('Phase')
ax.set_ylabel('ATLAS Relative Flux')
fig.savefig(output_dir+str(ra)+'_'+str(dec)+'_ATLAS_binned.png', dpi=300)
print('Saved '+output_dir+str(ra)+'_'+str(dec)+'_ATLAS_binned.png')

fig, ax  =plt.subplots(figsize=(10,4))
ax.errorbar(t, y, dy, fmt='.k', ms=2, elinewidth=1, capsize=1., alpha=0.8)
ax.set_xlabel('MBJD')
ax.set_ylabel('ATLAS Relative Flux')
fig.savefig(output_dir+str(ra)+'_'+str(dec)+'_ATLAS_raw.png', dpi=300)
print('Saved '+output_dir+str(ra)+'_'+str(dec)+'_ATLAS_raw.png')
