# Testing q factor range
# AM CVn TIC 1201247611

import numpy as np
import lc_utils as lcu
from Period_Finding import BLS, LS_Astropy

ticid = 1201247611
data_dir = "/home/echickle/data/s0057/s0057-lc/"

# Where to save BLS, LS plots 
output_dir = "/home/echickle/data/s0057/AMCVn_test/"
bls_dir    = output_dir+"bls"
ls_dir     = output_dir+"ls"
os.makedirs(bls_dir, exist_ok=True)
os.makedirs(ls_dir, exist_ok=True)

# Light curve extraction parameters
n_std = 3 
wind = 0.1

# BLS & LS period minimum & maximum
pmin = 400 / 60
pmax = 0.13

# ------------------------------------------------------------------------------

# Make text file 
ft = open(output_dir+"freq_table.txt", "w")
cols = []
ft.write("{}\t"*len(cols)+"\n")
ft.close()

# Get time and flux array
ticid_list = np.load(data_dir+'id-4-3.npy')
ind = np.nonzero(ticid_list == ticid)[0][0]
t = np.load(data_dir+'ts-4-3.npy')
y = np.load(data_dir+'lc-4-3.npy')[ind]

qrng = [[0.005, 0.05], [0.01, 0.05], [0.05,0.3]]

for i in range(len(qrng)):
    qmin, qmax = qrng[i]    
    t, y, flag = lcu.prep_lc(t, y, n_std=n_std, wind=wind)


    # Run BLS
    _, _, _, per, pow_best, bls_freq, bls_power, dur, epo = \
            BLS(t,y,dy,pmin=pmin,pmax=pmax,qmin=qmin,qmax=qmax,remove=True)
    

    _, _, _, per_ls, ls_power_best, ls_frq, ls_pow = LS_Astropy(t,y,dy,pmax=pmax)

    

