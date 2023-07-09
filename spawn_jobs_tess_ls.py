import os
import sys
import numpy as np
import multiprocessing
from multiprocessing import Pool
from run_ls import run_process


# # >> UZAY
# sector, cam, ccd = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
# output_dir = "/scratch/echickle/s%04d-WD/"%sector
# data_dir = "/scratch/data/tess/lcur/ffi/s%04d-lc/"%sector

# >> ENGAGING
N=int(sys.argv[1])
ccd_list = []
for sector in [56,57,58,59,60,61,62,63,64,65]: # 160 sublists
    for cam in [1,2,3,4]:
        for ccd in [1,2,3,4]:
            ccd_list.append([sector,cam,ccd])
# !! 
ccd_list = [[57,2,2], [58,1,4], [58,4,1], [59,2,1], [61,1,3], [61,1,4], 
            [61,3,2], [61,4,3], [62,1,2], [62,1,3], [62,1,4], [62,4,3],
            [63,1,1], [63,1,3], [63,1,4], [64,2,4], [65,4,3]]
sector, cam, ccd = ccd_list[N-1]
output_dir = "/pool001/echickle/tess/s%04d-WD/"%sector
data_dir = "/nobackup1c/users/echickle/TESS_Lightcurves/s%04d-lc/"%sector

# >> set up result directories
ls_dir = output_dir + "LS_cam{}-ccd{}/".format(cam,ccd)
gpu_dir = output_dir + "LS_results/"
os.makedirs(output_dir, exist_ok=True)
os.makedirs(ls_dir, exist_ok=True)
os.makedirs(gpu_dir, exist_ok=True)

# >> load ticid
suffix = "-{}-{}.npy".format(cam, ccd)
ticid = np.int64(np.load(data_dir+'id'+suffix))

# # >> remove completed 
# fnames_c = os.listdir(ls_dir)
# ticid_c = [str(f.split('_')[12][3:]) for f in fnames_c if f.split('.')[-1] == 'png']
# ticid_c = np.array(ticid_c)
# inter, comm1, comm2 = np.intersect1d(ticid, ticid_c, return_indices=True)
# ticid = np.delete(ticid, comm1) 

p = []
for i in range(len(ticid)):
    p.append( [sector, cam, ccd, ticid[i], data_dir, ls_dir] )

if __name__ == "__main__":
    # N_p=10 # >> uzay has 64 cores, 4 GPUs => 16 CPUs per GPU
    N_p = 16 # >> engaging
    multiprocessing.set_start_method('spawn')
    pool = Pool(processes=N_p)
    result = pool.map(run_process, p) # ticid, ra, dec, sig, snr, wid, period, period_min, q, phi0, epo, rp, nt, dphi
    lengths = [len(res) for res in result]
    print(np.unique(lengths))
    try:
        fmt='%s,%10.7f,%10.7f,%10.5f,%i,%10.5f,%10.8f,%10.5f'
        header='ticid, ra, dec, sig, wid, per, per_min, dphi'

        np.savetxt(gpu_dir+'LS-{}-{}-{}.result'.format(sector, cam, ccd), np.array(result),
                   fmt=fmt, header=header)
    except:
        print('np.savetxt failed!')
        file = open(gpu_dir+'LS-{}-{}-{}_unformatted.result'.format(sector, cam, ccd), 'w')
        for res in result:
            file.write(str(res)+'\n')
        file.close()

# f = '/pool001/echickle/tess/s0061-WD/s0061-gpu-res/GPU-61-3-2_unformatted.result'
# lines = f.readlines()
# lines = [line[1:-2].split(',') for line in lines]
# np.savetxt(f, np.array(lines).astype('float'), fmt='%s,%10.5f,%10.5f,%10.5f,%10.5f,%i,%10.8f,%10.5f,%10.5f,%10.5f,%10.5f,%10.5f,%i,%10.5f', header='ticid, ra, dec, sig, snr, wid, period, period_min, q, phi0, epo, rp, nt, dphi')
