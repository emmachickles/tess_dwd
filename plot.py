img_dir = '/data/submit/echickle/out/s41/'
out_dir = '/data/submit/echickle/out/s41-keck/'
data_dir = '/data/submit/echickle/data/sector_41_lc/'
# ticid_list = [1868421156, 1201294971, 1102498665, 1001328172, 1102501819,
#               1102431501, 1102539257, 1201319601, 1201272729, 1102517373,
#               335682765, 1102377238, 1102587443, 1201303785, 10015410257]

import numpy as np
import matplotlib.pyplot as plt
import lc_utils as lcu
# import fnmatch
import os

#ticid_list = [1201247611, 1883504789, 1400855795]
#ticid_list = [1201247611, 1883504789, 0]
ticid_list = [1201247611, 1400855795, 0]

# mag_list = [19, 19.5, 19]
# period_list = [49.7]
#period_list = [49.70946718807686, 226.71, 39.3427365936]
period_list = [49.70946718807686, 96.49, 39.3427365936]

#period_list = [49.70946718807686, 226.71, 96.49]
#t0_list = [20/1440, 150/1440, 150/1440]
t0_list = [20/1440, 4/1440, 150/1440]

# t0_list = [20/1440, 150/1440, 4/1440]

# # >> g
# mag = [18.12694, 18.50531, 19.06677]

# # >> bp
# mag = [18.05615, 18.58795, 19.05699]

# # >> rp
# mag = [17.83857, 18.29737, 18.63019]

# >> Pan-STARRS I
# mag_list = [19.029399871826172, 18.87820053100586, 19.27239990234375]
# mag_list = [19.029399871826172, 18.87820053100586, 15.6412]
mag_list = [19.029399871826172, 19.27239990234375, 15.6412]

ra_list = ['16:11:33.966', '19:36:19.946', '17:45:37.398']
dec_list = ['+63:08:31.67', '+54:09:20.74', '+60:19:24.21']

bls = os.listdir(img_dir)
qmin = 0.005
qmax = 0.2
bins = 120

fig, ax = plt.subplots(ncols=len(ticid_list), figsize=(10,3))

for i in range(len(ticid_list)):

    # -- prepare light curve ---------------------------------------------------

    if i ==2:
        data_dir2= '/data/submit/tess/echickle/s0056-lc/backup-2-2/'
        ticid = np.load(data_dir2+'id-2-2.npy')
        ind = np.nonzero(ticid == ticid_list[i])
        flux = np.load(data_dir2+'lc-2-2-ap0.7-bkg[1.3, 1.7].npy')
        time = np.load(data_dir2+'ts-2-2.npy')
        
    else:
        camccd = None
        while type(camccd) == type(None):
            for cam in [1,2,3,4]:
                for ccd in [1,2,3,4]:
                    suffix = '-'+str(cam)+'-'+str(ccd)+'.npy'
                    ticid = np.load(data_dir+'id'+suffix).astype('int')
                    ind = np.nonzero(ticid == ticid_list[i])
                    if len(ind[0]) > 0:
                        camccd = suffix

        suffix = camccd
        print(suffix)

        ticid = np.load(data_dir+'id'+suffix).astype('int')
        ind = np.nonzero(ticid == ticid_list[i])                    

        flux = np.load(data_dir+'lc'+suffix)
        time = np.load(data_dir+'ts'+suffix)
        
    t, y = lcu.prep_lc(time, flux[ind][0], n_std=3, wind=0.1)
    dy = np.ones(y.shape)*0.2
    
    # -- get best BLS fit ------------------------------------------------------

    # import cuvarbase.bls as bls
    # baseline = max(t) - min(t)
    # df = qmin / baseline
    # fmin = 1440 / (period_list[i] + 1)
    # fmax = 1440 / (period_list[i] - 1)
    # nf = int(np.ceil((fmax - fmin) / df))
    # freqs = fmin + df * np.arange(nf)
    
    # bls_power, sols = bls.eebls_gpu(t, y, dy, freqs, qmin=qmin, qmax=qmax)
    # q_best, phi0_best = sols[np.argmax(bls_power)]
    # f_best = freqs[np.argmax(bls_power)]

    # -- plot ------------------------------------------------------------------


    if i == 0:
        bins = 70
    else:
        bins=100
    folded_t, folded_y, folded_dy = lcu.make_phase_curve(t, y, period_list[i]/1440, bins=bins,
                                                         t0=t[0]+t0_list[i])
    shift = np.max(folded_t) - np.min(folded_t)
    ax[i].errorbar(folded_t*1440, folded_y, yerr=folded_dy*0.1,
                   fmt='.k', ms=1, elinewidth=1)
    # ax[i].errorbar((folded_t+shift)*1440, folded_y, yerr=folded_dy*0.1,
    #                fmt='.k', ms=1, elinewidth=1)
    ax[i].set_xlabel('Time [minutes]')
    if i ==2:
        ax[i].set_title('ZTF J2130+4420\nmax mag '+str(round(mag_list[i],2)))
    else:
        ax[i].set_title('TIC '+str(ticid_list[i])+'\nmax mag '+str(round(mag_list[i],2)))
    
    # fig, ax = plt.subplots()
    # lcu.plot_orbital_phase_curve(ax, t, y, dy, f_best, q_best, phi0_best)
    # fig.savefig('/home/submit/echickle/foo.png')
    
    # lcu.plot_phase_curve(t, y, period_list[i]/1440, out_dir, bins=100,
    #                      prefix='TIC%016d'%ticid_list[i]+'-')
    # lcu.plot_phase_curve(t, y, period_list[i]/1440, out_dir, bins=100,
    #                      prefix='TIC%016d'%ticid_list[i]+'-')

ax[0].set_ylabel('Normalized Flux')    
fig.tight_layout()
plt.savefig(out_dir+'phase_curve.png', dpi=300)
print('Saved '+out_dir+'phase_curve.png')
