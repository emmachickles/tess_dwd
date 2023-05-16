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

# Kevin\'s\ UCBs\ -\ UCBs.csv
# from KBB_Utils.KBB_Utils import LC_Tools
import LC_Tools
    
def run_lc_extraction(data_dir, out_dir, wd_cat, cam=None, ccd=None, mult_output=False,
                      tica=False):
    # >> load white dwarf catalog
    sources=np.loadtxt(wd_cat, usecols=(0,1,2), dtype='str')
    ticid_main=sources[:,0]
    catalog_main=SkyCoord(ra=sources[:,1].astype('float')*u.degree,
                          dec=sources[:,2].astype('float')*u.degree,
                          frame='icrs')
    if cam:
        cam_list = [int(cam)]
    else:
        cam_list = [1,2,3,4]

    if ccd:
        ccd_list = [int(ccd)]
    else:
        ccd_list = [1,2,3,4]

    os.makedirs(out_dir, exist_ok=True)

    # print(cam_list)
    # print(ccd_list)
    # with open('/home/submit/echickle/foo.txt', 'w') as f:
    #     f.write(str(cam_list))
    #     f.write('/n')
    #     f.write(str(ccd_list))

    for cam in cam_list:
        print('Running camera '+str(cam))
        for ccd in ccd_list:
            ccd_dir = 'cam{}-ccd{}/'.format(cam, ccd)
            
            p = [data_dir+ccd_dir+f for f in os.listdir(data_dir+ccd_dir)]
            run_ccd(p, catalog_main, ticid_main, cam, ccd, out_dir, 
                    mult_output=mult_output, tica=tica)
            print('Finished cam {} ccd {}!'.format(cam, ccd))

def run_target(data_dir, out_dir, cam, ccd, name, ra, dec):

    ticid_main=np.array([name])
    catalog_main=SkyCoord(ra=[ra*u.degree],
                          dec=[dec*u.degree],
                          frame='icrs')
    
    ccd_dir = 'cam{}-ccd{}/'.format(cam, ccd)
    suffix = '-'+name

    p = [data_dir+ccd_dir+f for f in os.listdir(data_dir+ccd_dir)]
    run_ccd(p, catalog_main, ticid_main, cam, ccd, out_dir, 
            mult_output=mult_output, suffix=suffix)
    
            
def run_UCBs():
    from astroquery.mast import Catalogs    
    # Kevin Burdge's Ultra-Compact Binaries

    fname = "/data/submit/echickle/Kevin\'s UCBs - UCBs.csv"

    with open(fname, 'r') as f:
        lines = f.readlines()[1:]
    id_main = []
    ra, dec = [], []
    for i in range(len(lines)):
        id_main.append(lines[i].split(',')[1])
        ra.append(float(lines[i].split(',')[2]))
        dec.append(float(lines[i].split(',')[3]))

        
        
    id_main, ra, dec = np.array(id_main), np.array(ra), np.array(dec)
    catalog_main=SkyCoord(ra=ra*u.degree,
                          dec=dec*u.degree,
                          frame='icrs')
        
    
    pass
            
def download_ccd(curl_file, data_dir, cam, ccd):
    ccd_dir = 'cam{}-ccd{}/'.format(cam, ccd)
    os.makedirs(data_dir+ccd_dir, exist_ok=True)
    with open(curl_file, 'r') as f:
        lines = f.readlines()[1:]
    for line in lines:
        ffi_cam = int(line.split(' ')[5].split('-')[2])
        ffi_ccd = int(line.split(' ')[5].split('-')[3])
        if ffi_cam == int(cam) and ffi_ccd == int(ccd):
            line = line.split(' ')
            line[5] = data_dir+ccd_dir+line[5]
            if not os.path.exists(line[5]):
                line = ' '.join(line)
                os.system(line)

def download_ccd_tica(sect_dir, sector, cam, ccd):
    for orbit in ['o1a', 'o1b', 'o2a', 'o2b']:
        curl_file=sect_dir+"hlsp_tica_tess_ffi_s%04d"%sector+\
            "-"+orbit+"-cam{}-ccd{}_tess_v01_ffis.sh".format(cam,ccd)
        with open(curl_file, 'r') as f:
            lines = f.readlines()[1:]
        for line in lines:
            line = line.split(' ')
            line[4] = sect_dir+line[4][1:-1]
            
            if not os.path.exists(line[4]):
                line = ' '.join(line)
                os.system(line)

def check_download(data_dir, cam, ccd):
    import subprocess
    ccd_dir = data_dir+'cam{}-ccd{}/'.format(cam, ccd)
    fnames = os.listdir(ccd_dir)
    fnames = np.array(fnames)
    sizes = []
    for i in range(len(fnames)):
        sizes.append(os.path.getsize(ccd_dir+fnames[i]))
    sizes = np.array(sizes)
    sbin, cnts = np.unique(sizes, return_counts=True)
    for i in range(len(sbin)):
        if sbin[i] < np.max(sizes) and cnts[i] < 0.1*len(fnames):
            inds = np.nonzero(sizes == sbin[i])
            for j in range(len(fnames[inds])):
                os.system('curl -C - -L -o '+ccd_dir+fnames[inds][j]+' https://mast.stsci.edu/api/v0.1/Download/file/?uri=mast:TESS/product/'+str(fnames[inds][j]))
    os.system('ls -lS {} | tail '.format(ccd_dir))

def check_download_tica(sect_dir, sector, cam, ccd):
    import subprocess

    # for orbit in ['o1a', 'o1b', 'o2a', 'o2b']:
    #     curl_file=sect_dir+"hlsp_tica_tess_ffi_s%04d"%sector+\
    #         "-"+orbit+"-cam{}-ccd{}_tess_v01_ffis.sh".format(cam,ccd)

    
        # curl_file=sect_dir+"hlsp_tica_tess_ffi_s%04d"%sector+\
        #     "-o1a-cam{}-ccd{}_tess_v01_ffis.sh".format(cam,ccd)
    
    ccd_dir = sect_dir+'s%04d/'%sector+'cam{}-ccd{}/'.format(cam, ccd)
    fnames = os.listdir(ccd_dir)
    fnames = np.array(fnames)
    sizes = []
    for i in range(len(fnames)):
        sizes.append(os.path.getsize(ccd_dir+fnames[i]))
    sizes = np.array(sizes)
    sbin, cnts = np.unique(sizes, return_counts=True)
    for i in range(len(sbin)):
        if sbin[i] < np.max(sizes) and cnts[i] < 0.001*len(fnames):
            inds = np.nonzero(sizes == sbin[i])
            for j in range(len(fnames[inds])):
                os.system('curl -f --create-dirs --output '+ccd_dir+fnames[inds][j]+\
                          ' https://mast.stsci.edu/api/v0.1/Download/file/?uri=mast:HLSP/tica/s%04d/'%sector+\
                          'cam{}-ccd{}/'.format(cam,ccd)+fnames[inds][j])
    os.system('ls -lS {} | tail '.format(ccd_dir))
                
def source_list(f,catalog,ticid, tica=False):

    hdu_list = fits.open(f)
    hd = hdu_list[0].header

    if tica:
        hd2=hdu_list[0].header
        t=hd['MJD-BEG']+100/86400
        # t = hd['STARTTJD'] # in TJD        
        dt = hdu_list[0].data    
    else:
        hd2=hdu_list[1].header
        t=hd['TSTART'] # BTJD
        dt = hdu_list[1].data        
        
    n_dty = dt.shape[0]
    n_dtx = dt.shape[1]
    w = wcs.WCS(hd2)
    print(w)
    
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
        t = hd['STARTTJD'] # in TJD
        cadence = hd['CADENCE']
        dt = hdu_list[0].data
    else:
        hd2=hdu_list[1].header # >> calibrated ffi
        cadence = hd['FFIINDEX']
        # TESS Barycentric Julian Day (BTJD), this is a Julian day minus 2457000.0 and corrected to the arrival times at the barycenter of the Solar System
        t=hd['TSTART'] # in BTJD

        dt = hdu_list[1].data # >> calibrated ffi 
    n_dty = dt.shape[0]
    n_dtx = dt.shape[1]
    w = wcs.WCS(hd2)

    
    #print(w)

    # image = hdu_list[0].data

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

    return t, cadence, phot_bkgsub


def run_ccd(p, catalog_main, ticid_main, cam, ccd, out_dir, mult_output=False,
            suffix='', tica=False):
    suffix = '-'+str(cam)+'-'+str(ccd)+suffix
    
    LC=[]
    ts=[]
    cadence_list=[]
    orbit_list = []
    
    trimmed_catalog, ticid, central_coord=source_list(p[0],catalog_main, ticid_main, tica=tica)

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
    failed_inds = []
    for f in p:
        try:
            t, cn, fluxes=process(f,aperture,background, central_coord, tica=tica)                
            ts.append(t)
            cadence_list.append(cn)
            LC.append(fluxes)
        except:
            print('Failed '+str(n_iter))
            failed_inds.append(n_iter)
        if n_iter // 10:
            print(n_iter)

        n_iter += 1

        # if n_iter // 10: # >> save intermediate products
        #     np.save(out_dir+'ts'+suffix+'.npy', np.array(ts))
        #     np.save(out_dir+'lc'+suffix+'.npy', np.array(LC).T)

    # -- save timeseries -------------------------------------------------------
    ts=np.array(ts)
    cadence_list = np.array(cadence_list)
    orbit_list = np.array(orbit_list)
    LC=np.array(LC)

    inds = np.argsort(ts)
    ts = ts[inds]
    cadence_list = cadence_list[inds]

    np.save(out_dir+'ts'+suffix+'.npy', ts)    
    np.save(out_dir+'cn'+suffix+'.npy', cadence_list)
   
    if mult_output:
        col = 0
        for i in range(len(N_ap_list)):
            for j in range(len(N_bkg_list)):
                suffix1 = '-ap'+str(N_ap_list[i])+'-in'+str(N_bkg_list[j][0])+\
                    '-out'+str(N_bkg_list[j][1])
                LC_p = LC[:,col,:].T
                LC_p = LC_p[:,inds]
                np.save(out_dir+'lc'+suffix+suffix1+'.npy', LC_p)
                col += 1
    else:
        LC=LC.T
        LC = LC[:,inds]
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

# ------------------------------------------------------------------------------

sector = 63

# >> file paths
# wd_cat    = '/home/echickle/data/WDs.txt'
wd_cat    = '/home/echickle/data/ZTF_Eclipses.txt'

sect_dir  = '/home/echickle/data/s%04d/'%sector
data_dir  = sect_dir+'s%04d/'%sector

# out_dir   = sect_dir+'s%04d-lc/'%sector
out_dir   = sect_dir+'s%04d-lc-ZTF/'%sector

curl_file = data_dir + 'tesscurl_sector_{}_ffic.sh'.format(sector)

# >> light curve parameters
N_ap  = 0.7
N_in  = 1.5
N_out = 2
mult_output = False # >> produce multiple light curves per source
N_ap_list   = [0.5, 0.7, 0.9, 1.1]
N_bkg_list  = [[1.3, 1.7], [1.8, 2.3], [1.8, 2.], [1.5, 2]]

if len(sys.argv) > 1:
    cam = sys.argv[1]
    ccd = sys.argv[2]

# -- RUN SETTINGS --------------------------------------------------------------
    
tica = False
cam = 2

for ccd in [1,2,3,4]:

    download_ccd(curl_file, data_dir, cam, ccd)
    check_download(data_dir, cam, ccd)
    try:
        run_lc_extraction(data_dir, out_dir, wd_cat, cam=cam, ccd=ccd,
                         mult_output=mult_output, tica=tica)
    except:
        pass
    os.system('rm -r '+data_dir+'cam{}-ccd{}'.format(cam,ccd))

    # download_ccd_tica(sect_dir, sector, cam, ccd)
    # check_download_tica(sect_dir, sector, cam, ccd)
    # run_lc_extraction(data_dir, out_dir, wd_cat, cam=cam, ccd=ccd,
    #                  mult_output=mult_output, tica=tica)
    # os.system('rm -r '+data_dir+'cam{}-ccd{}'.format(cam,ccd))
    
# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

# -- target run ----------------------------------------------------------------

# s, cam, ccd, name, ra, dec = 56, 2, 1, "ZTF J222827.07+494916.4", 337.1127, 49.82125
# s, cam, ccd, name, ra, dec = 56, 2, 2, "ZTF J213056.71+442046.5", 322.7362856, 44.34622882
# s, cam, ccd, name, ra, dec = 56, 3, 2, "ZTF J205515.98+465106.5", 313.8165708, 46.8517723
# s, cam, ccd, name, ra, dec = 57, 2, 2, "ZTF J224342.97+544206.1", 340.9290431, 52.70166019
# s, cam, ccd, name, ra, dec = 57, 2, 2, "ZTF J222827.07+494916.4", 337.1127, 49.82125
# s, cam, ccd, name, ra, dec = 57, 3, 2, "ZTF J235115.4+630527.7", 357.8141279, 63.0910333
# s, cam, ccd, name, ra, dec = 57, 2, 3, "ZTF J232020.43+375030.8", 350.0851861, 37.84185019

# data_dir   = '/data/submit/tess/echickle/s00{}/'.format(s)
# out_dir   = '/data/submit/tess/echickle/KBUCB-lc/'
# run_target(data_dir, out_dir, cam, ccd, name, ra, dec)


# ==============================================================================
