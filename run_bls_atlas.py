# -- inputs --------------------------------------------------------------------

pmin = 410 / 60 
pmax = 0.13 
qmin = 0.01
qmax = 0.15

data_dir = "/nobackup1c/users/echickle/ATLAS/"
output_dir = "/nobackup1c/users/echickle/out/"
bls_dir    = output_dir + 'bls'
ls_dir     = output_dir + 'ls'

fnames = ['1821925406145803136', '4590026743168925952', '6015541335204372864', '2674459092793700736', '6203449070679085312', '6448964626977687168', '6077174596944106752', '3946393595708149888', '5388880416727504640', '6884555017522418688', '4493071773470153088', '6673895637682011008', '1261495322912293248', '6806227939462649728', '4651860906336301696', '145831353928524672', '4405206054089274496', '5432401011098810368', '4906058133388497152', '5076475907343058304']

# ------------------------------------------------------------------------------

import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt

from Period_Finding import BLS, LS_Astropy
import lc_utils as lcu

import sys
import os

# ------------------------------------------------------------------------------

os.makedirs(output_dir, exist_ok=True)
os.makedirs(bls_dir, exist_ok=True)
os.makedirs(ls_dir, exist_ok=True)

# ------------------------------------------------------------------------------

for i in range(len(fnames)):
    if i % 50 == 0:
        print(str(i)+" / "+str(len(fnames)))

    t, y, dy = load_atlas_lc(data_dir + fnames[i])

    # -- compute BLS ----------------------------------------------------    
    t, y, dy, period, bls_power_best, freqs, power, q, phi0 = \
        BLS(t,y,dy,pmin=pmin,pmax=pmax,qmin=qmin,qmax=qmax,remove=True)
    prefix1 = 'ATLAS_'+prefix+'_'
    vet_plot(t, y, freqs, power, q, phi0, output_dir=bls_dir+prefix1,
             ticid=ticid, dy=dy)


    # -- compute LS ----------------------------------------------------

    _, _, _, ls_period, ls_power_best, ls_freqs, ls_power = \
        LS_Astropy(t,y,dy,pmax=pmax)

    lcu.vet_plot(t, y, ls_freqs, ls_power,output_dir=ls_dir+prefix1,
                 bins=100, save_npy=False, bls=False)

