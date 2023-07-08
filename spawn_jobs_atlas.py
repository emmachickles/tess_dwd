import pdb
import os
import sys
import numpy as np
import multiprocessing
from multiprocessing import Pool
from run_bls_atlas import run_process

# desc = "WDUCB"
# data_dir = "/home/echickle/data/atlas_"+desc+"/"
# data_dir = "/home/echickle/data/atlasforcedphotometryresults_"+desc+"/"
# output_dir = "/home/echickle/out/"
# bls_dir = output_dir + 'plot_'+desc+'/'

# data_dir = "/matchfiles/data2/ATLAS/"
# data_dir = "/home/echickle/data/atlasforcedphotometryresults_ZTF/"
# data_dir = "/home/echickle/data/atlasforcedphotometryresults/"

# data_dir = "/pool001/echickle/ATLAS/"
data_dir = "/pool001/echickle/sdB/"
output_dir = "/pool001/echickle/sdB_res/"
bls_dir = output_dir + 'bls_plot/'

os.makedirs(output_dir, exist_ok=True)
os.makedirs(bls_dir, exist_ok=True)

# N_p = 3 # >> hypernova.Caltech.edu
N_p=16 # >> engaging

p = [data_dir+f for f in os.listdir(data_dir)]

N=int(sys.argv[1])
# N_sub = 500 # 
N_sub = 25 # sdB
sub=int(np.ceil(len(p)/N_sub)) # 1393864 files, 500 sublists = 2780 files per job
if N == N_sub:
    p=p[(N-1)*sub:]
else:    
    p=p[(N-1)*sub:N*sub]

# p = [data_dir+f for f in os.listdir(data_dir)]
# p = p[:20]
# p = [data_dir+str(6079447764213678592)]
# p = []

# >> Kevin's UCBs
# gid = [282679289838317184, 2931430078489330944, 2134541781964072320, 1916879054217244544,
#        1334928852673292544, 3402486457838073728, 3223663226718611328, 1805551543400732544,
#        4488756396492138880, 1869111286948848128]

# >> TESS eclipsing systems
# gid = [3070712053964938624, 3043990932113003648, 578539413395848704, 580790014913812608,
#        3077510098136276480, 3084580919978268672, 5509633937648970624, 5579580606800791808,
#        4654758944077464192, 5665509162794290944, 5670594640294856576, 5673444097693647616,
#        5621623354466444800, 5473641183995972608, 5468670738602933504, 5463274408550748160,
#        5462557110355830912, 5302515565071142400, 5242787486412627072, 5210507882302442368]
# for i in range(len(gid)):
#     p.append(data_dir+str(gid[i]))

for i in range(len(p)):
    p[i] = [p[i], data_dir, bls_dir]


if __name__ == '__main__':
    multiprocessing.set_start_method('spawn')
    pool = Pool(processes=N_p)
    result = pool.map(run_process, p) # gaiaid, ra, dec, sig, snr, wid, period, period_min, q, phi0, epo, rp, nt, dphi

    # np.savetxt(output_dir+'GPU_'+desc+'.result',np.array(result),
    #            fmt='%s', header='gaiaid ra dec sig snr wid period period_min q phi0 epo rp nt dphi')
    # print(output_dir+'GPU_'+desc+'.result')

    np.savetxt(output_dir+'GPU_'+str(N)+'.result',np.array(result),
               fmt='%s', header='gaiaid ra dec sig snr wid period period_min q phi0 epo rp nt dphi')
    print(output_dir+'GPU_'+str(N)+'.result')
