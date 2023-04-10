import numpy as np
import matplotlib.pyplot as plt
import lc_utils as lcu
import os

tess_dir = "/scratch/data/tess/lcur/ffi/"
ticid  = [2040677137]
sector = [57]
cam    = [2]
ccd    = [3]
period = [0.038365738]

for i in range(len(ticid)):
    data_dir = tess_dir + "s%04d-lc/"%sector[i]
    suffix = '-{}-{}.npy'.format(cam[i], ccd[i])

    t = np.load(data_dir+'ts'+suffix)
    ticid_list = np.load(data_dir+'id'+suffix)
    ind = np.nonzero(ticid_list == ticid[i])[0][0]
    y = np.load(data_dir+'lc'+suffix)[ind]
    t, y, flag =lcu.prep_lc(t, y)

    folded_t, folded_y, folded_dy = lcu.bin_timeseries(t%period[i], y, bins=100)
    fig, ax = plt.subplots()
    lcu.plot_phase_curve(ax, folded_t, folded_y, folded_dy, period[i],
                     ylabel="TESS Relative Flux")
    fig.savefig("/scratch/echickle/MGAB/"+str(ticid[i])+".png")
    print("/scratch/echickle/MGAB/"+str(ticid[i])+".png")

    
