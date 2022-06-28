#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 26 11:14:04 2022

@author: kburdge
"""

import numpy as np
import matplotlib.pyplot as plt
import pathlib

from astropy import wcs
from astropy.io import fits
from photutils import aperture_photometry, CircularAperture, SkyCircularAnnulus
from photutils import SkyCircularAperture

from astropy import units as u
from astropy.coordinates import SkyCoord

# >> edited by etc 220603
import pdb
import sys
import os
sys.path.append(os.getcwd())
from KBB_Utils.KBB_Utils import LC_Tools
import matplotlib.pyplot as plt

# sector = sys.argv[1]
# CCD = sys.argv[2]
# quad = sys.argv[3]

# p = [str(pp) for pp in pathlib.Path('Sector'+str(sector)+'/').glob('*-'+str(int(CCD))+'-'+str(int(quad))+'-*.fits')]

data_dir = '/scratch/echickle/data/sector_30_ffic/'
p = [data_dir+f for f in os.listdir(data_dir)]
sources=np.loadtxt('tess_dwd/WDs.txt',usecols=(1,2))
# <<


catalog_main=SkyCoord(ra=sources[:,0]*u.degree, dec=sources[:,1]*u.degree, frame='icrs')

def source_list(f,catalog):

    hdu_list = fits.open(f)
    hd = hdu_list[0].header
    #print(hd)
    hd2=hdu_list[1].header
    t=hd['TSTART']+300/86400.0
    dt = hdu_list[1].data
    #print(dt)
    n_dty = dt.shape[0]
    n_dtx = dt.shape[1]
    w = wcs.WCS(hd2)
    
    central_coord=SkyCoord.from_pixel(1024,1024,w)
    
    print(central_coord)
    
    #idxc, idxcatalog, d2d, d3d = catalog.search_around_sky(central_coord, 10*u.deg)

    
    d2d = central_coord.separation(catalog)
    catalogmsk = d2d < np.sqrt(2)*6*u.deg
    idxcatalog = np.where(catalogmsk)[0]
    
    
    trimmed_catalog=catalog[catalogmsk]
    sky_aperture = SkyCircularAperture(trimmed_catalog, r=21 * u.arcsec)
    
    
    pix_aperture = sky_aperture.to_pixel(w)
    
    print(pix_aperture)
    mask=np.abs(pix_aperture.positions[:,0]>0) & np.abs(pix_aperture.positions[:,1]>0) & np.abs(pix_aperture.positions[:,0]<2048) & np.abs(pix_aperture.positions[:,1]<2048)
    trimmed_catalog2=trimmed_catalog[mask]
    


    
    #print(trimmed_catalog)


    return trimmed_catalog2

def process(f,sky_aperture,background):

    hdu_list = fits.open(f)
    hd = hdu_list[0].header
    #print(hd)
    hd2=hdu_list[1].header
    t=hd['TSTART']+300/86400.0
    dt = hdu_list[1].data
    #print(dt)
    n_dty = dt.shape[0]
    n_dtx = dt.shape[1]
    w = wcs.WCS(hd2)
    
    #print(w)

    image = hdu_list[0].data

    image = np.array(dt)
    #image -= np.nanmedian(dt)


    pix_aperture = sky_aperture.to_pixel(w)
    phot_table = aperture_photometry(image, pix_aperture)  
    
    background_pix_aperture = background.to_pixel(w)
    background_phot_table = aperture_photometry(image, background_pix_aperture)  
    
    norm=background_pix_aperture.area/pix_aperture.area
    print(norm)
    
    for col in phot_table.colnames:  
        phot_table[col].info.format = '%.8g'  # for consistent table out


    return t, phot_table['aperture_sum'].value-background_phot_table['aperture_sum'].value/norm
LC=[]
ts=[]
trimmed_catalog=source_list(p[0],catalog_main)

N_ap=0.7
aperture = SkyCircularAperture(trimmed_catalog, r=N_ap*21 * u.arcsec)
N_in=1.5
N_out=2
background = SkyCircularAnnulus(trimmed_catalog, N_in*21*u.arcsec,N_out*21*u.arcsec)
# p=p[:2] # >> etc 220603
for f in p:
    #coord= SkyCoord(ra=242.891526160*u.degree, dec=63.142131440*u.degree, frame='icrs')
    
    try:
        t, fluxes=process(f,aperture,background)
        ts.append(t)
        LC.append(fluxes)
    except:
        pass

ts=np.array(ts)    
LC=np.array(LC)
print(LC)
LC=LC.T

# >> edited by etc 220603
np.save('/scratch/echickle/dwd/ts.npy', ts)
np.save('/scratch/echickle/dwd/lc.npy', LC)

ticid_sources=np.loadtxt('tess_dwd/WDs.txt',usecols=(0))
ticid = []
for i in range(len(trimmed_catalog)):
    idx, d2d, d3d = trimmed_catalog[i].match_to_catalog_sky(catalog_main)
    ticid.append(ticid_sources[idx])
np.save('/scratch/echickle/dwd/id.npy', ticid)

lc = np.load('/scratch/echickle/dwd/lc.npy').T
ts = np.load('/scratch/echickle/dwd/ts.npy')

os.makedirs('/scratch/echickle/dwd/raw_lcs')
for i in range(len(lc)):
    plt.figure()
    plt.plot(ts, lc[i], '.k')
    plt.xlabel('Time')
    plt.ylabel('Flux')
    plt.savefig('/scratch/echickle/dwd/raw_lcs/lc_'+str(i)+'.png')
    plt.close()
    print(i)
# <<
