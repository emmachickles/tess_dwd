import numpy as np
import matplotlib.pyplot as plt
import lc_utils as lcu
import os

tess_dir = "/scratch/data/tess/lcur/ffi/"
ticid  = [808963820]
sector = [63]
cam    = [3]
ccd    = [2]
period = [0.038365738/2]

for i in range(len(ticid)):
    data_dir = tess_dir + "s%04d-lc/"%sector[i]
    suffix = '-{}-{}.npy'.format(cam[i], ccd[i])

    t = np.load(data_dir+'ts'+suffix)
    ticid_list = np.load(data_dir+'id'+suffix)
    ind = np.nonzero(ticid_list == ticid[i])[0][0]
    y = np.load(data_dir+'lc'+suffix)[ind]

    inds = np.nonzero( t> 59870 )
    t, y = t[inds],y[inds]

    t, y, flag =lcu.prep_lc(t, y)

    folded_t, folded_y, folded_dy = lcu.bin_timeseries(t%period[i], y, bins=100)
    fig, ax = plt.subplots()
    lcu.plot_phase_curve(ax, folded_t, folded_y, folded_dy, period[i],
                     ylabel="TESS Relative Flux")
    fig.savefig("/scratch/echickle/KB_UCB/"+str(ticid[i])+".png")
    print("/scratch/echickle/KB_UCB/"+str(ticid[i])+".png")

    
