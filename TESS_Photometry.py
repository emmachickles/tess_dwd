#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 26 11:14:04 2022

@author: kburdge
"""
import pdb
import sys
import os
sys.path.append('/home/submit/echickle/work/')

import numpy as np
import matplotlib.pyplot as plt
import pathlib

from astropy    import wcs
from astropy    import units as u
from astropy.io import fits
from photutils  import aperture_photometry, CircularAperture, SkyCircularAnnulus
from photutils  import SkyCircularAperture
from astropy.coordinates import SkyCoord

from KBB_Utils.KBB_Utils import LC_Tools

wd_cat    = '/data/submit/echickle/WDs.txt'
data_dir  = '/data/submit/tess/echickle/s0057/'
out_dir   = '/data/submit/tess/echickle/s0057-lc/'
curl_file = '/data/submit/echickle/data/tesscurl_sector_41_ffic.sh'

# >> light curve parameters
N_ap  = 0.7
N_in  = 1.5
N_out = 2

mult_output = False # >> produce multiple light curves per source
N_ap_list   = [0.5, 0.7, 0.9, 1.1]
N_bkg_list  = [[1.3, 1.7], [1.8, 2.3], [1.8, 2.], [1.5, 2]]

cam = sys.argv[1]

def run_sector(data_dir, out_dir, wd_cat, curl_file, cam=None, mult_output=False):
    # >> load white dwarf catalog
    sources=np.loadtxt(wd_cat, usecols=(0,1,2))
    ticid_main=sources[:,0].astype('int')
    catalog_main=SkyCoord(ra=sources[:,1]*u.degree, dec=sources[:,2]*u.degree,
                          frame='icrs')
    if cam:
        cam_list = [int(cam)]
    else:
        cam_list = [1,2,3,4]

    for cam in cam_list:
        print('Running camera '+str(cam))
        for ccd in [1,2,3,4]:
            ccd_dir = 'cam{}-ccd{}/'.format(cam, ccd)
            
            p = [data_dir+ccd_dir+f for f in os.listdir(data_dir+ccd_dir)]
            run_ccd(p, catalog_main, ticid_main, cam, ccd, out_dir, wd_cat,
                    mult_output=mult_output)
            print('Finished cam {} ccd {}!'.format(cam, ccd))

def download_sector(data_dir, out_dir, wd_cat, curl_file, cam=None):
    if cam:
        cam_list = [int(cam)]
    else:
        cam_list = [1,2,3,4]

    for cam in cam_list:
        print('Running camera '+str(cam))
        for ccd in [1,2,3,4]: 
            download_ccd(curl_file, data_dir, cam, ccd)
            print('Downloaded cam {} ccd {}!'.format(cam, ccd))
            os.system('rm '+data_dir+ccd_dir+'*.fits')
            os.rmdir(data_dir+ccd_dir)
            
def download_ccd(curl_file, data_dir, cam, ccd):
    ccd_dir = 'cam{}-ccd{}/'.format(cam, ccd)
    os.makedirs(data_dir+ccd_dir)            
    with open(curl_file, 'r') as f:
        lines = f.readlines()[1:]
    for line in lines:
        ffi_cam = int(line.split(' ')[5].split('-')[2])
        ffi_ccd = int(line.split(' ')[5].split('-')[3])
        if ffi_cam == cam and ffi_ccd == ccd:            
            line = line.split(' ')
            line[5] = data_dir+ccd_dir+line[5]
            line = ' '.join(line)
            os.system(line)

def source_list(f,catalog,ticid, tica=True):

    hdu_list = fits.open(f)
    hd = hdu_list[0].header

    if tica:
        hd2=hdu_list[0].header
        t=hd['MJD-BEG']+100/86400        
        dt = hdu_list[0].data    
    else:
        hd2=hdu_list[0].header
        t=hd['TSTART']+300/86400.0
        dt = hdu_list[1].data        
        
    n_dty = dt.shape[0]
    n_dtx = dt.shape[1]
    w = wcs.WCS(hd2)
    print(w)
    
    central_coord=SkyCoord.from_pixel(1024,1024,w)
    t = BJDConvert(t,central_coord.ra.deg, central_coord.dec.deg,
                   date_format='mjd').value
    
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
    trimmed_ticid=ticid[catalogmsk][mask]

    return trimmed_catalog2, trimmed_ticid, central_coord

def BJDConvert(times, RA, Dec, date_format='mjd', telescope='Palomar'):
    '''Function for converting a series of timestamps to Barycentric Julian
    Date format in Barycentric Dynamical time'''
    import numpy as np
    from astropy.time import Time
    from astropy.coordinates import EarthLocation
    from astropy.coordinates import SkyCoord  # High-level coordinates
    from astropy.coordinates import ICRS, Galactic, FK4, FK5, BarycentricTrueEcliptic  # Low-level frames
    from astropy.coordinates import Angle, Latitude, Longitude  # Angles
    import astropy.units as u
    
    t = Time(times,format=date_format,scale='utc')
    t2=t.tcb
    c = SkyCoord(RA,Dec, unit="deg")
    d=c.transform_to(BarycentricTrueEcliptic)
    Observatory=EarthLocation.of_site(telescope)
    delta=t2.light_travel_time(c,kind='barycentric',location=Observatory)
    BJD_TCB=t2+delta
    return BJD_TCB

def get_flux(sky_aperture, background, wcs, image):

    from photutils import ApertureStats
    
    pix_aperture = sky_aperture.to_pixel(wcs=wcs)
    pix_annulus = background.to_pixel(wcs=wcs)

    sky_stats = ApertureStats(image, pix_aperture)
    bkg_stats = ApertureStats(image, pix_annulus)

    sky_area = sky_stats.sum_aper_area.value
    bkg_area = bkg_stats.sum_aper_area.value

    phot_bkgsub = sky_stats.sum - bkg_stats.median * sky_area    
    # phot_bkgsub = sky_stats.sum - bkg_stats.mean * sky_area
    # norm = bkg_area / sky_area
    # phot_bkgsub = sky_stats.sum - bkg_stats.sum / norm # background-subtracted flux
    return phot_bkgsub
    
def process(f,sky_aperture,background,central_coord, tica=True):    
    hdu_list = fits.open(f)
    hd = hdu_list[0].header
    #print(hd)

    if tica:
        hd2=hdu_list[0].header    
        t=hd['MJD-BEG']+100/86400
        dt = hdu_list[0].data
    else:
        hd2=hdu_list[1].header
        t=hd['TSTART']+300/86400.0
        dt = hdu_list[1].data        
    #print(dt)
    n_dty = dt.shape[0]
    n_dtx = dt.shape[1]
    w = wcs.WCS(hd2)

    t = BJDConvert(t,central_coord.ra.deg, central_coord.dec.deg, date_format='mjd').value    
    
    #print(w)

    image = hdu_list[0].data

    image = np.array(dt)
    #image -= np.nanmedian(dt)

    if type(sky_aperture) == type([]): # >> produce multiple light curves for each source
        phot = []
        for i in range(len(sky_aperture)):
            for j in range(len(background)):                

                phot = get_flux(sky_aperture[i], background[j], w, image)
                phot_bkgsub.append(phot)
        phot_bkgsub = np.array(photbkgsub)
        
    else: # >> produce single light curve for each source
        phot_bkgsub = get_flux(sky_aperture, background, w, image)
        
    return t, phot_bkgsub
                
def run_ccd(p, catalog_main, ticid_main, cam, ccd, out_dir, wd_cat, mult_output=False):
    suffix = '-'+str(cam)+'-'+str(ccd)
    
    LC=[]
    ts=[]
    
    trimmed_catalog, ticid, central_coord=source_list(p[0],catalog_main, ticid_main)

    if mult_output:
        aperture = []
        background = []
        for i in range(len(N_ap_list)):
            aperture.append(SkyCircularAperture(trimmed_catalog, r=N_ap_list[i]*21*u.arcsec))
        for j in range(len(N_bkg_list)):
            background.append(SkyCircularAnnulus(trimmed_catalog, N_bkg_list[j][0]*21*u.arcsec,
                                                 N_bkg_list[j][1]*21*u.arcsec))
    else:
        aperture = SkyCircularAperture(trimmed_catalog, r=N_ap*21 * u.arcsec)
        background = SkyCircularAnnulus(trimmed_catalog,
                                        N_in*21*u.arcsec,N_out*21*u.arcsec)    
    n_iter = 0
    for f in p:
        t, fluxes=process(f,aperture,background, central_coord)
        ts.append(t)
        LC.append(fluxes)
        n_iter += 1
        if n_iter // 10:
            print(n_iter)

        # if n_iter // 10: # >> save intermediate products
        #     np.save(out_dir+'ts'+suffix+'.npy', np.array(ts))
        #     np.save(out_dir+'lc'+suffix+'.npy', np.array(LC).T)

    # -- save timeseries -------------------------------------------------------
    ts=np.array(ts)    
    LC=np.array(LC)

    np.save(out_dir+'ts'+suffix+'.npy', ts)    
    
    if mult_output:
        col = 0
        for i in range(len(N_ap_list)):
            for j in range(len(N_bkg_list)):
                suffix1 = '-ap'+str(N_ap_list[i])+'-in'+str(N_bkg_list[j][0])+\
                    '-out'+str(N_bkg_list[j][1])
                np.save(out_dir+'lc'+suffix+suffix1+'.npy', LC[:,col,:].T)
                col += 1
    else:
        LC=LC.T
        np.save(out_dir+'lc'+suffix+'.npy', LC)

    # >> save ticid
    np.save(out_dir+'id'+suffix+'.npy', ticid)

    # >> save coordinates
    ra = []
    dec = []
    for i in range(len(trimmed_catalog)):
        ra.append(trimmed_catalog[i].ra.degree)
        dec.append(trimmed_catalog[i].dec.degree)
    co = np.array([ra, dec]).T
    np.save(out_dir+'co'+suffix+'.npy', co)        

run_sector(data_dir, out_dir, wd_cat, curl_file, cam=cam, mult_output=mult_output)

# ==============================================================================

# -- getting ra dec ------------------------------------------------------------
# sources=np.loadtxt(wd_cat, usecols=(0,1,2))
# ticid_main=sources[:,0].astype('int')
# catalog_main=SkyCoord(ra=sources[:,1]*u.degree, dec=sources[:,2]*u.degree,
#                       frame='icrs')
# for cam in [3]:
#     for ccd in [4]:
#         ticid_ccd = np.load(out_dir + 'id-{}-{}.npy'.format(cam, ccd))
#         ra, dec = [], []
#         for i in range(len(ticid_ccd)):
#             ind = np.nonzero(ticid_main == ticid_ccd[i])[0][0]
#             ra.append(catalog_main[ind].ra.degree)
#             dec.append(catalog_main[ind].dec.degree)
#         co = np.array([ra, dec]).T
#         np.save(out_dir+'co-{}-{}.npy'.format(cam, ccd), co)
        
# -- plotting ------------------------------------------------------------------
# for i in range(len(LC)):
#     plt.figure()
#     plt.plot(ts, LC[i], '.k')    
#     plt.title('N_ap: {}, N_in: {}, N_out: {}'.format(N_ap, N_in, N_out))
#     plt.xlabel('Time')
#     plt.ylabel('Flux')
#     plt.savefig(out_dir+'lc_'+str(i)+'.png', dpi=300)
#     plt.close()
#     print(i)

                # pix_aperture = sky_aperture[i].to_pixel(w)
                # phot_table = aperture_photometry(image, pix_aperture)  
                # background_pix_aperture = background[j].to_pixel(w)
                # background_phot_table = aperture_photometry(image, background_pix_aperture)  

                # norm=background_pix_aperture.area/pix_aperture.area
                # print(norm)

                # for col in phot_table.colnames:  
                #     phot_table[col].info.format = '%.8g'  # for consistent table out

                # phot.append(phot_table['aperture_sum'].value-background_phot_table['aperture_sum'].value/norm)

                        # pix_aperture = sky_aperture.to_pixel(w)
        # phot_table = aperture_photometry(image, pix_aperture)  

        # background_pix_aperture = background.to_pixel(w)
        # background_phot_table = aperture_photometry(image, background_pix_aperture)  

        # norm=background_pix_aperture.area/pix_aperture.area
        # print(norm)

        # for col in phot_table.colnames:  
        #     phot_table[col].info.format = '%.8g'  # for consistent table out

        # return t, phot_table['aperture_sum'].value-background_phot_table['aperture_sum'].value/norm

    # trimmed_catalog, ticid=source_list('/data/submit/tess/echickle/s0056/cam1ccd1/hlsp_tica_tess_ffi_s0056-o1-00690247-cam1-ccd1_tess_v01_img.fits',catalog_main, ticid_main)
        
    # flag = True
    # i = 0
    # while flag:
    #     try:
    #         trimmed_catalog, ticid, central_coord=source_list(p[i],catalog_main, ticid_main)
    #     except:
    #         i+=1
    #     else:
    #         flag = False

    # !!
    # if cam == 2 and ccd == 2:
    #     co = SkyCoord(322.7362855032135,44.346236059843456, unit='deg')
    #     trimmed_catalog = SkyCoord(list(trimmed_catalog) + [co])
    #     ticid = np.append(ticid, 0)

        # trimmed_catalog = SkyCoord([trimmed_catalog[0], co])
        # ticid= [ticid[0], 0000]
    
    # ticid_sources=np.loadtxt(wd_cat,usecols=(0))
    # ticid = []
    # for i in range(len(trimmed_catalog)):
    #     idx, d2d, d3d = trimmed_catalog[i].match_to_catalog_sky(catalog_main)
    #     ticid.append(ticid_sources[idx])
