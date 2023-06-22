import os
import multiprocessing
from multiprocessing import Pool
import lc_utils as lcu
import numpy as np

data_dir = '/data/ATLAS/'
output_dir = '/home/echickle/data/EOF_IQR1010/'
gpu_dir = output_dir + 'GPU_res/'
vet_dir = output_dir + 'vet_res/'
os.makedirs(vet_dir, exist_ok=True)

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
        pow, wid = cat[i][1], cat[i][2]

        if os.path.exists(data_dir+str(gid)):
            t, y, dy, ra, dec = lcu.load_atlas_lc(data_dir+str(gid),
                                                      pos_iqr=10, neg_iqr=10)
            snr, _, _, _, _, nt, dphi = lcu.calc_snr(t, y, per, q, phi0)
            res.append([gid, ra, dec, pow, snr, wid, per, q, phi0, nt, dphi])
    np.savetxt(vet_dir+f, np.array(res), fmt='%s')    
        

if __name__ == '__main__':
    N_p = 64
    multiprocessing.set_start_method('spawn')
    pool = Pool(processes=N_p)
    result = pool.map(run_process, p) 
    


