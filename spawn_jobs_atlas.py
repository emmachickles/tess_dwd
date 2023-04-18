import pdb
import os
import sys
import numpy as np
import multiprocessing
from multiprocessing import Pool
from run_bls_atlas import run_process

# data_dir = "/home/echickle/ATLAS_TEST/"
data_dir = "/pool001/echickle/ATLAS/"
p = [data_dir+f for f in os.listdir(data_dir)]

N=int(sys.argv[1])
N_sub = 500
sub=int(np.ceil(len(p)/N_sub)) # 1393864 files, 500 sublists = 2780 files per job
if N<N_sub-1:
    p=p[N*sub:(N+1)*sub]
else:
    p=p[N*sub:]

# p = p[:9] # !!
# p.insert(0, '/matchfiles/data2/ATLAS/3084580919978268672')

if __name__ == '__main__':
    N_p=16
    multiprocessing.set_start_method('spawn')
    pool = Pool(processes=N_p)
    result = pool.map(run_process, p) # gaiaid, sig, snr, wid, period, period_min, q, phi0, dur, epo
    np.savetxt('GPU'+str(N)+'.result',np.array(result),fmt='%i,%10.5f,%10.5f,%i,%10.5f,%10.5f,%10.5f,%10.5f,%10.5f,%10.5f')
