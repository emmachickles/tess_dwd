from lc_utils import get_tess_lc, normalize_lc, bin_timeseries, plot_phase_curve, prep_lc, rm_qflag
import numpy as np
import matplotlib.pyplot as plt
import pdb

out_dir = "/home/echickle/out/vet/"
tess_dir = "/home/echickle/data/"
qflag_dir = "/home/echickle/data/QLPqflags/"

bins=150
ztf=False
n_std=10
figsize=(5,2)

# >> 60 minute eclipsers
ra_list = [215.68014, 167.49633, 158.683, 123.67249]
dec_list = [-44.07947, -40.0399, -57.698669, -64.44789]
per_list = [0.04719063, 0.04393835, 0.04360338, 0.04412959]

ra_list = [285.3559028]
dec_list = [53.1581301]
per_list = [0.0281957]

bkg = 'w'
c= 'k'
if bkg == 'k':
    import matplotlib as mpl
    plt.rcParams['font.family'] = 'sans-serif'
    mpl.rcParams['text.color'] = 'white'
    mpl.rcParams['axes.labelcolor'] = 'white'
    mpl.rcParams['xtick.color'] = 'white'
    mpl.rcParams['ytick.color'] = 'white'
    mpl.rcParams['axes.edgecolor'] = 'white'
    

for ra, dec, per in zip(ra_list, dec_list, per_list):
    # ra, dec = 121.174886, -2.262535
    suffix='ra_{}_dec_{}'.format(ra,dec)
    ticid, sector, cam, ccd = get_tess_lc(tess_dir, ra=ra, dec=dec, ztf=ztf)
    # ticid = 938779482
    # sector = 64
    # cam = 2
    # ccd = 1
    if ztf:
        ccd_dir = tess_dir+'s00{}/s00{}-lc-ZTF/'.format(sector, sector)
    else:
        ccd_dir = tess_dir+'s00{}/s00{}-lc/'.format(sector, sector)
        
    f_suffix = '-{}-{}.npy'.format(cam, ccd)        
    t = np.load(ccd_dir+'ts'+f_suffix)
    cn = np.load(ccd_dir+'cn'+f_suffix)    
    ticid_list = np.load(ccd_dir+'id'+f_suffix)
    ticid, ticid_list = np.int64(ticid), np.int64(ticid_list)
    ind = np.nonzero(ticid_list == ticid)[0][0]
    y = np.load(ccd_dir+'lc'+f_suffix)[ind]

    fig00, ax00=plt.subplots(figsize=figsize)
    if bkg == 'k':
        fig00.patch.set_facecolor('black')
        ax00.set_facecolor('black')
    ax00.plot(t, y, '.'+c, ms=1)
    ax00.set_xlabel('Time [TJD]')
    ax00.set_ylabel('Background-subtracted flux')
    plt.tight_layout()
    fig00.savefig(out_dir+suffix+'_raw_TESS.png', dpi=300)
    print('Saved '+out_dir+suffix+'_raw_TESS.png')

    t, y, cn = rm_qflag(t, y, cn, qflag_dir, sector, cam, ccd)    
    t, y, flag = prep_lc(t, y, n_std=10)
    
    y, _ = normalize_lc(y)

    fig00, ax00=plt.subplots(figsize=figsize)    
    if bkg == 'k':
        fig00.patch.set_facecolor('black')
        ax00.set_facecolor('black')
        
    ax00.plot(t, y, '.'+c, ms=1)
    ax00.set_xlabel('Time [TJD]')
    ax00.set_ylabel('Relative flux')
    plt.tight_layout()
    fig00.savefig(out_dir+suffix+'_norm_TESS.png', dpi=300)
    print('Saved '+out_dir+suffix+'_norm_TESS.png')    

    folded_t, folded_y, folded_dy = bin_timeseries(t%per, y, bins)

    
    fig00, ax00=plt.subplots(figsize=figsize)

    if bkg == 'k':
        fig00.patch.set_facecolor('black')
        ax00.set_facecolor('black')
    
    folded_t, folded_y, folded_dy = bin_timeseries(t%per, y, bins)
    plot_phase_curve(ax00, folded_t, folded_y, folded_dy,
                     ylabel="TESS Relative Flux")
    # plot_phase_curve(ax00, folded_t, folded_y, folded_dy, period=per,
    #                  ylabel="TESS Relative Flux")    
    fig00.savefig(out_dir+suffix+'_binned_TESS.png', dpi=300)
    print('Saved '+out_dir+suffix+'_binned_TESS.png')

    


# # ------------------------------------------------------------------------------

# lc_dir = "/scratch/echickle/grnd_lc/"

# lc_f = lc_dir + "GPU-56-1-1.result" # !!
# sector, cam = 56, 1 # !!

# out_dir = "/scratch/echickle/LDSS_230421_ATLAS/"
# # out_dir = "/scratch/echickle/s%04d/"%sector \
# #           +"s%04d-"%sector+str(cam)+"-crossmatch/"
# tess_dir = "/scratch/data/tess/lcur/ffi/s%04d-lc/"%sector

# gaia_tab = "/scratch/echickle/100pc_clean.fits"
# bls=True

# # ------------------------------------------------------------------------------

# import lc_utils as lcu
# import numpy as np
# import sys
# import os
# import pdb

# # ------------------------------------------------------------------------------

# lc_tab = np.loadtxt(lc_f, delimiter=",", dtype="str", skiprows=1)
# if len(lc_tab.shape) == 1:
#     lc_tab = np.expand_dims(lc_tab, 0)
# ticid = np.unique(lc_tab[:,0])
# os.makedirs(out_dir, exist_ok=True)

# pmin = 400 / 60 
# pmax = 10
# qmin = 0.01
# qmax = 0.15

# # >> remove completed
# # fnames = os.listdir(out_dir)
# # ticid_out = [int(f.split('_')[13][3:]) for f in fnames if f.split('_')[0] == 'TESS']
# # ticid_out = np.array(ticid_out)
# # inter, comm1, comm2 = np.intersect1d(ticid, ticid_out, return_indices=True)
# # ticid = np.delete(ticid, comm1) 

# # !!
# # ticid = ['767706310']

# for i in range(len(ticid)):
#     print(str(i)+'/'+str(len(ticid)))
#     lc_info = lc_tab[lc_tab[:,0] == ticid[i]]
#     fname_tess = lc_info[lc_info[:,1] == 'TESS'][0][2]
#     fname_atlas = lc_info[lc_info[:,1] == 'ATLAS'][:,2]
#     if len(fname_atlas) > 0:
#         fname_atlas = lc_dir+fname_atlas[0]
#     else:
#         fname_atlas = None
#     fnames_ztf = list(lc_info[lc_info[:,1] == 'ZTF'][:,2])
#     for j in range(len(fnames_ztf)):
#         fnames_ztf[j] = lc_dir+fnames_ztf[j]

#     fname_tess = fname_tess.split('_')    
#     if bls:
#         ticid = int(fname_tess[12][3:])
#         cam = fname_tess[15]
#         ccd = fname_tess[17]
#         per = float(fname_tess[7]) / 1440        
#         ra = float(fname_tess[19])
#         dec = float(fname_tess[21])    
#         suffix = '_'.join(fname_tess[12:22])
#     else:
#         ticid = int(fname_tess[4][3:])
#         cam = fname_tess[6]
#         ccd = fname_tess[8]
#         per = float(fname_tess[3]) / 1440        
#         ra = float(fname_tess[10])
#         dec = float(fname_tess[12])    
#         suffix = '_'.join(fname_tess[4:13])

#     lcu.make_panel_plot(fname_atlas,fnames_ztf,tess_dir,ticid,cam,ccd,per,ra,dec,
#                         gaia_tab,out_dir,suffix,bls=bls,pmin=pmin,pmax=pmax,qmin=qmin,qmax=qmax)
