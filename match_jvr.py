data_dir = "/scratch/data/tess/lcur/ffi/s0063-lc-ZTF/"
# data_dir = "/home/echickle/data/s0062/s0062-lc-ZTF/"
import os
import numpy as np
from astroquery.mast import Catalogs

for cam in [1]:
    for ccd in [4]:
        suffix = '-{}-{}.npy'.format(cam, ccd)
        # suffix = suffix_list[i]

        old_fname = data_dir+'id'+suffix
        new_fname = data_dir+'zid'+suffix
        os.system('cp '+old_fname+' '+new_fname)

        co = np.load(data_dir+'co'+suffix)
        ticid = []
        for i in range(len(co)):
            cat = Catalogs.query_region('{} {}'.format(co[i][0], co[i][1]),
                                        radius=0.01, catalog='Tic')
            ticid.append( int(cat['ID'][0]) )

        np.save(data_dir+'id'+suffix, np.array(ticid, dtype='int'))
