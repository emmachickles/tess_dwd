import os
import sys
import numpy as np
import multiprocessing
from multiprocessing import Pool
from run_bls_tess import run_process

sector, cam, ccd = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])

# output_dir = "/scratch/echickle/s%04d/"%sector
output_dir = "/scratch/echickle/s%04d-ZTF/"%sector

# data_dir = "/scratch/data/tess/lcur/ffi/s%04d-lc/"%sector
data_dir = "/scratch/data/tess/lcur/ffi/s%04d-lc-ZTF/"%sector

bls_dir = output_dir + "s%04d"%sector + "-bls-{}-{}/".format(cam,ccd)
gpu_dir = output_dir + "s%04d"%sector + "-gpu-res/"
os.makedirs(output_dir, exist_ok=True)
os.makedirs(gpu_dir, exist_ok=True)

suffix = "-{}-{}.npy".format(cam, ccd)
ticid = np.load(data_dir+'id'+suffix).astype('int')

# >> remove completed 
# fnames_c = os.listdir(bls_dir)
# ticid_c = [str(f.split('_')[12][3:]) for f in fnames_c if f.split('.')[-1] == 'png']
# ticid_c = np.array(ticid_c)
# inter, comm1, comm2 = np.intersect1d(ticid, ticid_c, return_indices=True)
# ticid = np.delete(ticid, comm1) 

p = []
for i in range(len(ticid)):
    p.append( [sector, cam, ccd, ticid[i], data_dir, bls_dir] )

if __name__ == "__main__":
    N_p=16 # >> uzay has 64 cores, 4 GPUs => 16 CPUs per GPU
    multiprocessing.set_start_method('spawn')
    pool = Pool(processes=N_p)
    result = pool.map(run_process, p) # ticid, ra, dec, sig, snr, wid, period, period_min, q, phi0, dur, epo, rp, nt, dphi
    np.savetxt(gpu_dir+'GPU-{}-{}-{}.result'.format(sector, cam, ccd), np.array(result),
               fmt='%s,%10.5f,%10.5f,%10.5f,%10.5f,%i,%10.5f,%10.5f,%10.5f,%10.5f,%10.5f,%10.5f,%10.5f,%i,%10.5f',
               header='ticid, ra, dec, sig, snr, wid, period, period_min, q, phi0, dur, epo, rp, nt, dphi')
