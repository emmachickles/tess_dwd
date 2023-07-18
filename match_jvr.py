import sys
sector = int(sys.argv[1])

# data_dir = "/scratch/data/tess/lcur/ffi/s0056-lc-ZTF/" # !!
data_dir = "/home/echickle/data/s00{}/s00{}-lc-ZTF/".format(sector,sector)
import os
import numpy as np
from astroquery.mast import Catalogs

suffix_list = [f[2:] for f in os.listdir(data_dir) if 'co-' in f]
print(suffix_list)

for suffix in suffix_list:
    old_fname = data_dir+'id'+suffix
    new_fname = data_dir+'zid'+suffix # !!
    os.system('cp '+old_fname+' '+new_fname)

    co = np.load(data_dir+'co'+suffix)
    ticid = []
    for i in range(len(co)):
        if i % 100 == 0:
            print(i)
        cat = Catalogs.query_region('{} {}'.format(co[i][0], co[i][1]),
                                    radius=0.01, catalog='Tic')
        if len(cat) > 0:
            ticid.append( int(cat['ID'][0]) )
        else: 
            ticid.append(None)

    np.save(data_dir+'id'+suffix, np.array(ticid, dtype='int'))
