# -- inputs -------------------------------------------------------------------

pmin = 410 / 60 
pmax = 0.13 
qmin = 0.01
qmax = 0.15

data_dir = "/pool001/echickle/ATLAS/"
output_dir = "/nobackup1c/users/echickle/out/"
bls_dir    = output_dir + 'bls'

import sys
fname = sys.argv[1]

# ------------------------------------------------------------------------------

import matplotlib as mpl
mpl.use('Agg')

from Period_Finding import BLS
import lc_utils as lcu

# ------------------------------------------------------------------------------

os.makedirs(output_dir, exist_ok=True)
os.makedirs(bls_dir, exist_ok=True)

# ------------------------------------------------------------------------------

t, y, dy = lcu.load_atlas_lc(data_dir + fnames[i])

# -- compute BLS ---------------------------------------------------------------
t, y, dy, period, bls_power_best, freqs, power, q, phi0 = \
    BLS(t,y,dy,pmin=pmin,pmax=pmax,qmin=qmin,qmax=qmax,remove=True)
prefix1 = 'ATLAS_'+str(fnames[i])
lcu.vet_plot(t, y, freqs, power, q, phi0, output_dir=bls_dir+prefix1, dy=dy)
