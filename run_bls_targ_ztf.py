import lc_utils as lcu
import matplotlib.pyplot as plt
from Period_Finding import BLS
import numpy as np

qmin = 0.01
qmax = 0.15
output_dir = "/home/echickle/out/"

data_dir = '/home/echickle/'
# fnames_ztf = [data_dir+'251.78027_-24.91717_g.lc',
#               data_dir+'251.78027_-24.91717_r.lc']
# suffix = "251.78027_-24.91717"

ra = 275.26688
dec = 4.69748

fnames_ztf = lcu.get_ztf_lc(ra, dec)
suffix = str(ra)+"_"+str(dec)

period = 58.37516607/1440.


# -- load data -----------------------------------------------------------------

t, y, dy = lcu.load_ztf_lc(fnames_ztf, clip=False)

# -- run bls -------------------------------------------------------------------

# t, y, dy, period, bls_power_best, freqs, power, q, phi0 = \
#     BLS(t,y,dy,pmin=2,pmax=10,qmin=qmin,qmax=qmax,remove=False)
# lcu.vet_plot(t, y, freqs, power, q, phi0, output_dir=output_dir,
#              dy=dy, suffix=suffix)

# -- plot phase curves ---------------------------------------------------------


fig, ax = plt.subplots()
folded_t, folded_y, folded_dy = lcu.bin_timeseries(t%period, y, 100, dy=dy)
lcu.plot_phase_curve(ax, folded_t, folded_y, folded_dy, period=period,
                 ylabel="ZTF Flux")
fig.savefig(output_dir+suffix+'_binned.png', dpi=300)
print(output_dir+suffix+'_binned.png')

fig, ax=plt.subplots()
lcu.plot_phase_curve(ax, t%period, y, dy, period=period,
                     ylabel="ZTF Flux", alpha=0.6)
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

