import lc_utils as lcu
from Period_Finding import BLS


qmin = 0.01
qmax = 0.15
output_dir = "/scratch/echickle/"
fnames_ztf = ['260.16736_43.78046_g.lc',
              '260.16736_43.78046_i.lc',
              '260.16736_43.78046_r.lc']
suffix="_260.16736_43.78046"

t, y, dy = lcu.load_ztf_lc(fnames_ztf)
t, y, dy, period, bls_power_best, freqs, power, q, phi0 = \
    BLS(t,y,dy,pmin=2,pmax=10,qmin=qmin,qmax=qmax,remove=False)
lcu.vet_plot(t, y, freqs, power, q, phi0, output_dir=output_dir,
             dy=dy, plot_threshold=0, suffix=suffix)

