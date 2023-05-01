import lc_utils as lcu
from Period_Finding import BLS

qmin = 0.01
qmax = 0.15
output_dir = "/home/echickle/out/"
data_dir = "/matchfiles/data2/ATLAS/"
fname_atlas = data_dir+'5369709061413007360'
suffix="_5369709061413007360"

t, y, dy = lcu.load_atlas_lc(fname_atlas)
t, y, dy, period, bls_power_best, freqs, power, q, phi0 = \
    BLS(t,y,dy,pmin=2,pmax=10,qmin=qmin,qmax=qmax,remove=False)
lcu.vet_plot(t, y, freqs, power, q, phi0, output_dir=output_dir,
             dy=dy, suffix=suffix)
