import os
import multiprocessing
from multiprocessing import Pool
import lc_utils as lcu
import numpy as np

data_dir = '/data/ATLAS/'
output_dir = '/home/echickle/data/EOF_IQR1010/'
gpu_dir = output_dir + 'GPU_res/'
vet_dir = output_dir + 'vet_plot/'
os.makedirs(vet_dir, exist_ok=True)


pos_iqr = 10
neg_iqr = 10
skiprows = 0
objid_type = 'GAIAID'

wd_main = "/data/GaiaEDR3_WD_main.fits"
rp_ext = "/data/GaiaEDR3_WD_RPM_ext.fits"

pow_threshold=25
snr_threshold=1.
per_threshold=210
wid_threshold=5


# >> metrics = ['pow', 'snr', 'wid', 'nt', 'dphi']

p = os.listdir(gpu_dir)
# p = p[:4]

def run_process(f):
    # >> cols = ['gid', 'pow', 'wid', 'per', 'q', 'phi0']
    cat = np.loadtxt(gpu_dir+f, dtype='float', usecols=(0,1,3,4,6,7), delimiter=',')

    res = []
    for i in range(len(cat)):
        gid = np.int64(cat[i][0])
        per, q, phi0 = cat[i][3], cat[i][4], cat[i][5]
        sig, wid = cat[i][1], cat[i][2]

        f_lc = data_dir+str(gid)
        if os.path.exists(f_lc): 
            t, y, dy, ra, dec = lcu.load_atlas_lc(f_lc, pos_iqr=pos_iqr, neg_iqr=neg_iqr,
                                                  skiprows=skiprows)

            suffix = '_GID_'+str(gid)+'_ra_{}_dec_{}_'.format(ra, dec)    
            lcu.vet_plot(t, y, q=q, phi0=phi0, dy=dy, output_dir=vet_dir, suffix=suffix,
                         objid=gid, objid_type=objid_type,
                         ra=ra, dec=dec, sig=sig, wid=wid, period=per,
                         wd_main=wd_main, rp_ext=rp_ext,
                         pow_threshold=pow_threshold, snr_threshold=snr_threshold,
                         per_threshold=per_threshold, wid_threshold=wid_threshold)
        

if __name__ == '__main__':
    N_p = 64
    multiprocessing.set_start_method('spawn')
    pool = Pool(processes=N_p)
    result = pool.map(run_process, p)
    # run_process(p[0])
    


