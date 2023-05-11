import lc_utils as lcu
from Period_Finding import BLS

qmin = 0.01
qmax = 0.15
output_dir = "/home/echickle/out/"
# data_dir = "/matchfiles/data2/ATLAS/"
data_dir = '/data/ATLAS/'
gid = 580790014913812608
fname_atlas = data_dir + str(gid)
suffix="_"+str(gid)

t, y, dy, ra, dec = lcu.load_atlas_lc(fname_atlas, clip=False)
print(ra)
print(dec)
t, y, dy, period, bls_power_best, freqs, power, q, phi0 = \
    BLS(t,y,dy,pmin=2,pmax=10,qmin=qmin,qmax=qmax,remove=False)
lcu.vet_plot(t, y, freqs, power, q, phi0, output_dir=output_dir,
             dy=dy, suffix=suffix)
