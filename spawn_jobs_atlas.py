import os
import multiprocessing
from multiprocessing import Pool

data_dir = "/pool001/echickle/ATLAS/"
p = os.listdir(data_dir)

N=int(sys.argv[1])
quarter=int(np.ceil(len(p)/4.0))
print(quarter)
if N<3:
    p=p[N*quarter:(N+1)*quarter]
else:
    p=p[N*quarter:]

N_p=8
multiprocessing.set_start_method('spawn')
pool = Pool(processes=N_p)
pool.map(run_bls_atlas, p)
