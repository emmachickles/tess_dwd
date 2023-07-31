import os 
import pdb
import numpy as np
from vet_utils import *

data_dir = '/scratch/echickle/tess/BLS_results/'
out_dir = '/scratch/echickle/tess/metric_cut/'
plot_dir = '/scratch/echickle/tess/plots/'
sector_list = [56,57,58,59,60,61,62,63,64,65]
per_max = 60
nt_min = 50
dphi_max = 0.03
pow_min = 25
wid_min = 6
# plots already thresholded by pow > 25, wid > 6

# Set up output directory
os.makedirs(out_dir, exist_ok=True)

# Load plot names
plot_fnames = []
plot_ticid = []
for sector in sector_list:
    for cam in [1,2,3,4]:
        for ccd in [1,2,3,4]:
            ccd_dir=plot_dir+'s00{}-WD/BLS_cam{}-ccd{}/'.format(sector,cam,ccd)
            fnames = os.listdir(ccd_dir)
            ticid = [np.int64(z.split('_')[8][3:]) for z in fnames]
            fnames = [ccd_dir+f for f in fnames]
            plot_fnames.extend(fnames)
            plot_ticid.extend(ticid)
plot_fnames = np.array(plot_fnames)
plot_ticid = np.int64(np.array(plot_ticid))

# Load Gaia white dwarf results 
# cols : ticid, ra, dec, power, snr, wid, per, nt, dphi
result_list = append_result_file(data_dir, sector_list)

# Apply metric cuts
good_idx = np.nonzero( (result_list[:,6]<per_max/1440.) * (result_list[:,7]>nt_min) * (result_list[:,8]<dphi_max) * \
                       (result_list[:,3]>pow_min) * (result_list[:,5]>wid_min))
result_list = result_list[good_idx]

# Remove known soures
for i, ticid in enumerate(result_list[:,0]):
    if i % 100 == 0:
        print(i)
    ind = np.nonzero(plot_ticid == np.int64(ticid))[0]
    if len(ind) > 0:
        fname = str(plot_fnames[ind][0])
        os.system('cp '+fname+' '+out_dir)
pdb.set_trace()
