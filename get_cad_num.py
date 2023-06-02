# Get cadence numbers given BTJD observation start times
# Set up to run on hypernova.mit.edu
import os
from astropy.io import fits
import numpy as np
import pdb

sector = 61
data_dir = '/home/echickle/data/s%04d/'%sector
curl_f = data_dir + 's%04d/tesscurl_sector_%d_ffic.sh'%(sector,sector)
# 245700

f = open(curl_f, 'r')
lines = f.readlines()[1:]
dates = [l.split(' ')[5][4:17] for l in lines]

# >> download ffi from first orbit
ind = np.argsort(dates)[0]
line = lines[ind]
os.system(line)

# >> get time and cadence number from ffi
ffic = line.split(' ')[5]
hdul = fits.open(ffic)
tstart_o1 = hdul[0].header['TSTART'] 
ffiindex_o1 = hdul[0].header['FFIINDEX']

# >> download ffi from second orbit
ind = np.argsort(dates)[-1]
line = lines[ind]
os.system(line)

# >> get time and cadence number from ffi
ffic = line.split(' ')[5]
hdul = fits.open(ffic)
tstart_o2 = hdul[0].header['TSTART'] 
ffiindex_o2 = hdul[0].header['FFIINDEX']

# >> match
cadence = 200. / (24*60*60)
for cam in [1,2,3,4]:
    for ccd in [1,2,3,4]:
        ts=np.load(data_dir+'s%04d-lc/ts-%d-%d.npy'%(sector,cam,ccd))
        pdb.set_trace()
        cn=[]
        
        # >> get orbit gap
        ind = np.argmax(np.diff(np.sort(ts)))
        tgap = np.sort(ts)[ind] # >> second orbit has t > tgap

        for t in ts:
            if t <= tgap: # >> orbit 1
                cn.append(ffiindex_o1 + int( (t - tstart_o1) / cadence))
            else:
                cn.append(ffiindex_o2 + int( (t - tstart_o2) / cadence))

        np.save(data_dir+'s%04d-lc/cn-%d-%d.npy'%(sector,cam,ccd), cn)
