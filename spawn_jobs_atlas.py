import pdb
import os
import sys
import numpy as np
import multiprocessing
from multiprocessing import Pool
from run_bls_atlas import run_process

data_dir = "/home/echickle/ATLAS_TEST/"
# data_dir = "/pool001/echickle/ATLAS/"
p = [data_dir+f for f in os.listdir(data_dir)]

N=int(sys.argv[1])
quarter=int(np.ceil(len(p)/4.0))
print(quarter)
if N<3:
    p=p[N*quarter:(N+1)*quarter]
else:
    p=p[N*quarter:]

p = p[:17] # !!

if __name__ == '__main__':
    N_p=4
    multiprocessing.set_start_method('spawn')
    pool = Pool(processes=N_p)
    pool.map(run_process, p)
