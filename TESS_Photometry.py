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
import pandas as pd
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
                      tica=False, save_dir=None):
    if wd_cat.split('/')[-1] == 'Gaia_blue.csv':
        sources=pd.read_csv(wd_cat)
        ticid_main=sources['source_id'].to_numpy()
        catalog_main=SkyCoord(ra=sources['ra'].to_numpy()*u.degree,
                              dec=sources['dec'].to_numpy()*u.degree, frame='icrs')
    else:
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
                    mult_output=mult_output, tica=tica, save_dir=save_dir)
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
        # t=hd['MJD-BEG']+100/86400
        t = hd['STARTTJD'] # in TJD        
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
    if len(trimmed_catalog) == 0:
        return None, None, None

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

    sky_area = np.pi * pix_aperture.r**2
    bkg_area = np.pi * (pix_annulus.r_out**2 - pix_annulus.r_in**2)
    
    sky_stats = ApertureStats(image, pix_aperture)
    bkg_stats = ApertureStats(image, pix_annulus)

    # sky_area = sky_stats.sum_aper_area.value
    # bkg_area = bkg_stats.sum_aper_area.value

    # phot_bkgsub = sky_stats.sum - bkg_stats.median * sky_area    
    # phot_bkgsub = sky_stats.sum - bkg_stats.mean * sky_area
    norm = bkg_area / sky_area
    phot_bkgsub = sky_stats.sum - bkg_stats.sum / norm # background-subtracted flux
    # if np.count_nonzero(np.isnan(phot_bkgsub)) > 0:
    #     pdb.set_trace()
    return phot_bkgsub


def visualize_flux(sky_aperture, background, wcs, image, phot_bkgsub, save_dir, suffix=''):
    from photutils import SkyCircularAperture
    from astropy.visualization import ZScaleInterval
    from astropy.coordinates import SkyCoord
    from astroquery.gaia import Gaia
    from photutils import ApertureStats
    
    pix_aperture = sky_aperture.to_pixel(wcs=wcs)
    pix_annulus = background.to_pixel(wcs=wcs)

    # sky_stats = ApertureStats(image, pix_aperture)
    # bkg_stats = ApertureStats(image, pix_annulus)

    # sky_area = sky_stats.sum_aper_area.value
    # bkg_area = bkg_stats.sum_aper_area.value

    # phot_bkgsub = sky_stats.sum - bkg_stats.median * sky_area    

    # Plot the image
    fig, ax = plt.subplots()
    interval = ZScaleInterval()
    vmin, vmax = interval.get_limits(image)
    ax.imshow(image, cmap='gray', origin='lower', vmin=vmin, vmax=vmax)

    # Overlay circular apertures on the image
    radius = pix_aperture.r
    for i, position in enumerate(pix_aperture.positions):
        aperture = plt.Circle(position, radius, edgecolor='cyan', facecolor='none', lw=2.)
        ax.add_patch(aperture)
        # ax.text(position[0], position[1], 'Flux {:5f}'.format(phot_bkgsub[i]), color='cyan', fontsize=8, ha='left', va='bottom')
        ax.text(position[0], position[1], 'RA {} Dec {}\nFlux {:5f}'.format(sky_aperture.positions[i].ra.degree, sky_aperture.positions[i].dec.degree, phot_bkgsub[i]), color='cyan', fontsize=5, ha='left', va='bottom')

    # r_in = pix_annulus.r_in
    # r_out = pix_annulus.r_out
    # for position in pix_annulus.positions:
    #     aperture = plt.Circle(position, r_in, edgecolor='blue', facecolor='none', lw=0.5)
    #     ax.add_patch(aperture)
    #     aperture = plt.Circle(position, r_out, edgecolor='blue', facecolor='none', lw=0.5)        
    #     ax.add_patch(aperture)

    ax.set(title='Apertures')

    # Save the figure to the specified directory
    if not os.path.exists(save_dir):
        os.makedirs(save_dir)
    plt.savefig(os.path.join(save_dir, 'flux_visualization'+suffix+'.png'), dpi=300)
    plt.close()

    # Calculate the number of rows and columns for the subplot grid
    num_positions = len(pix_aperture.positions)
    num_columns = 4  # Specify the desired number of columns
    num_rows = (num_positions + num_columns - 1) // num_columns

    # Plot the image
    fig, axs = plt.subplots(num_rows, num_columns, figsize=(3.5*num_columns, 3.5 * num_rows))

    for i, position in enumerate(pix_aperture.positions):
        row = i // num_columns
        col = i % num_columns
        ax = axs[row, col]

        ax.imshow(image, cmap='gray', origin='lower', vmin=vmin, vmax=vmax)

        # Overlay circular apertures on the image
        aperture = plt.Circle(position, radius, edgecolor='red', facecolor='none', lw=2)
        ax.add_patch(aperture)

        r_in = pix_annulus.r_in
        r_out = pix_annulus.r_out
        annulus_inner = plt.Circle(position, r_in, edgecolor='blue', facecolor='none', lw=2)
        ax.add_patch(annulus_inner)
        annulus_outer = plt.Circle(position, r_out, edgecolor='blue', facecolor='none', lw=2)
        ax.add_patch(annulus_outer)

        ax.set(title='RA {} Dec {}\nFlux {:5f}'.format(sky_aperture.positions[i].ra.degree, sky_aperture.positions[i].dec.degree, phot_bkgsub[i]))

        # Zoom in on the image within the subplot
        x_min = position[0] - 1.15 * r_out
        x_max = position[0] + 1.15 * r_out
        y_min = position[1] - 1.15 * r_out
        y_max = position[1] + 1.15 * r_out
        ax.set_xlim(x_min, x_max)
        ax.set_ylim(y_min, y_max)
        
        ax.axis('off')


        # Add pixel values as text in the middle of each pixel
        for y in range(int(y_min), int(y_max)+1):
            for x in range(int(x_min), int(x_max+1)):
                value = image[y, x]
                ax.text(x, y, '%.2f'%value, color='fuchsia', ha='center', va='center', fontsize=7)
                
        # Convert pixel coordinates to celestial coordinates (RA and Dec)
        coord = SkyCoord.from_pixel(position[0], position[1], wcs)        
        
        # Query Gaia for sources in the field of view
        width = 1.5 * background.r_out
        height = 1.5 * background.r_out
        gaia_sources = Gaia.query_object(coord, width=width, height=height)
        gaia_sources.sort('phot_g_mean_mag')

        pdb.set_trace()
        # Plot Gaia sources as text on the subplot
        for source in gaia_sources[:5]:
            ra = source['ra']
            dec = source['dec']
            gid = source['source_id']
            co = SkyCoord(ra=ra, dec=dec, unit='deg') 
            
            pix = co.to_pixel(wcs=wcs)
            ax.plot(pix[0], pix[1], 'o', color='cyan')
            ax.text(pix[0], pix[1], 'GID {:d}'.format(gid), color='cyan', fontsize=10, ha='left', va='bottom')

        
    # Remove any empty subplots
    if num_positions < num_rows * num_columns:
        for j in range(num_positions, num_rows * num_columns):
            row = j // num_columns
            col = j % num_columns
            fig.delaxes(axs[row, col])

    # Adjust spacing between subplots
    # fig.tight_layout()
    # plt.subplots_adjust(wspace=0.2, hspace=0.3)        
    plt.savefig(os.path.join(save_dir, 'flux_visualization_subplot'+suffix+'.png'))
    print(os.path.join(save_dir, 'flux_visualization_subplot'+suffix+'.png'))
    pdb.set_trace()
    plt.close()

# def barycentric_correction_test(f, save_dir):
#     hdu_list = fits.open(f)
#     tstart = hdu_list[1].header['TSTART']
#     barycorr = hdu_list[1].header['BARYCORR']
#     tsat = tstart - barycorr
#     btc_col = hdu_list[1].header['BTC_PIX1'] # reference col for barycentric time correction
#     btc_row = hdu_list[1].header['BTC_PIX2'] # reference row for barycentric time correction
    
def process(f,sky_aperture,background,central_coord, tica=True, save_dir=None):    
    hdu_list = fits.open(f)
    hd = hdu_list[0].header

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
    print(w)
    image = np.array(dt)
    #image -= np.nanmedian(dt)

    if type(sky_aperture) == type([]): # >> produce multiple light curves for each source
        phot_bkgsub = []
        for i in range(len(sky_aperture)):
            for j in range(len(background)):                

                phot = get_flux(sky_aperture[i], background[j], w, image)
                phot_bkgsub.append(phot)

                if save_dir is not None:
                    visualize_flux(sky_aperture[i], background[i], w, image, phot, save_dir,
                                   suffix='_{:.5f}_{:.5f}_r_{:.2f}_rin_{:.2f}_rout_{:.2f}_{:.5f}'.format(central_coord.ra.degree, central_coord.dec.degree, sky_aperture[i].r.value, background[i].r_in.value, background[i].r_out.value, t))
                
        phot_bkgsub = np.array(phot_bkgsub)
        
    else: # >> produce single light curve for each source
        phot_bkgsub = get_flux(sky_aperture, background, w, image)
        if save_dir is not None:
            visualize_flux(sky_aperture, background, w, image, phot_bkgsub, save_dir, suffix='_{:.5f}_{:.5f}_{:.5f}'.format(central_coord.ra.degree, central_coord.dec.degree, t))
        

    return t, cadence, phot_bkgsub


def run_ccd(p, catalog_main, ticid_main, cam, ccd, out_dir, mult_output=False,
            suffix='', tica=False, save_dir=None, max_sources=15000):
    suffix = '-'+str(cam)+'-'+str(ccd)+suffix
    
    LC=[]
    ts=[]
    cadence_list=[]
    orbit_list = []
    
    trimmed_catalog, ticid, central_coord=source_list(p[0],catalog_main, ticid_main, tica=tica)
    
    if trimmed_catalog is None:
        return

    if len(trimmed_catalog) == 0:
        return

    if len(trimmed_catalog) > max_sources:
        trimmed_catalog_list = []
        n_split = len(trimmed_catalog) // max_sources
        for i in range( n_split ):
            trimmed_catalog_list.append(trimmed_catalog[i*max_sources:(i+1)*max_sources])
        if len(trimmed_catalog) % max_sources > 0:
            trimmed_catalog_list.append( trimmed_catalog[n_split*max_sources:])
    else:
        trimmed_catalog_list = [trimmed_catalog]
    
    for n, trimmed_catalog in enumerate(trimmed_catalog_list):
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
        for i, f in enumerate(p): 
            # t, cn, fluxes=process(f,aperture,background, central_coord, tica=tica, save_dir=save_dir) 
            # ts.append(t)
            # cadence_list.append(cn)
            # LC.append(fluxes)

            try:
                t, cn, fluxes=process(f,aperture,background, central_coord, tica=tica, save_dir=save_dir)  
                ts.append(t)
                cadence_list.append(cn)
                LC.append(fluxes)
            except:
                # pdb.set_trace()
                print('Failed '+str(i))
                failed_inds.append(i)

            print(i)

            # if n_iter // 10: # >> save intermediate products
            #     np.save(out_dir+'ts'+suffix+'.npy', np.array(ts))
            #     np.save(out_dir+'lc'+suffix+'.npy', np.array(LC).T)


        print('Number failed: '+str(len(failed_inds)))
        # -- save timeseries -------------------------------------------------------
        ts=np.array(ts)
        cadence_list = np.array(cadence_list)
        orbit_list = np.array(orbit_list)
        LC=np.array(LC)

        inds = np.argsort(ts)
        ts = ts[inds]
        cadence_list = cadence_list[inds]

        if len(trimmed_catalog_list) > 1:
            np.save(out_dir+'ts'+suffix+'_'+str(n)+'.npy', ts)    
            np.save(out_dir+'cn'+suffix+'_'+str(n)+'.npy', cadence_list)
        else:            
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
                    if len(trimmed_catalog_list) > 1:
                        np.save(out_dir+'lc'+suffix+suffix1+'_'+str(n)+'.npy', LC_p)
                    else:
                        np.save(out_dir+'lc'+suffix+suffix1+'.npy', LC_p)
                    col += 1
        else:
            LC=LC.T
            LC = LC[:,inds]
            if len(trimmed_catalog_list)> 1:
                np.save(out_dir+'lc'+suffix+'_'+str(n)+'.npy', LC)
            else:
                np.save(out_dir+'lc'+suffix+'.npy', LC)

        # >> save ticid
        if len(trimmed_catalog_list) > 1:
            np.save(out_dir+'id'+suffix+'_'+str(n)+'.npy', ticid)
        else:
            np.save(out_dir+'id'+suffix+'.npy', ticid)

        # >> save coordinates
        ra = []
        dec = []
        for i in range(len(trimmed_catalog)):
            ra.append(trimmed_catalog[i].ra.degree)
            dec.append(trimmed_catalog[i].dec.degree)
        co = np.array([ra, dec]).T
        if len(trimmed_catalog_list) > 1:
            np.save(out_dir+'co'+suffix+'_'+str(n)+'.npy', co)
        else:
            np.save(out_dir+'co'+suffix+'.npy', co)        

# ------------------------------------------------------------------------------

sector = 56

sect_dir  = '/home/echickle/data/s%04d/'%sector
data_dir  = sect_dir+'s%04d/'%sector

# >> file paths
# wd_cat    = '/home/echickle/data/WDs.txt'
# wd_cat    = '/home/echickle/data/ZTF_Eclipses.txt'

# out_dir   = sect_dir+'s%04d-lc/'%sector
# out_dir   = sect_dir+'s%04d-lc-ZTF/'%sector
# os.makedirs(out_dir, exist_ok=True)

# plot_dir  = sect_dir+'plot/'
# os.makedirs(plot_dir, exist_ok=True)

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
cam = 1

# ccd = 2
# run_lc_extraction(data_dir, out_dir, wd_cat, cam=cam, ccd=ccd,
#                   mult_output=mult_output, tica=tica)#, save_dir=plot_dir)

for ccd in [1,2,3,4]: # !!

    if tica:
        download_ccd_tica(sect_dir, sector, cam, ccd)    
        check_download_tica(sect_dir, sector, cam, ccd)
    else:
        download_ccd(curl_file, data_dir, cam, ccd)
        check_download(data_dir, cam, ccd)    

    # wd_cat    = '/home/echickle/data/WDs.txt'
    # out_dir   = sect_dir+'s%04d-lc/'%sector
    # os.makedirs(out_dir, exist_ok=True)
    # run_lc_extraction(data_dir, out_dir, wd_cat, cam=cam, ccd=ccd,
    #                  mult_output=mult_output, tica=tica)

    # wd_cat    = '/home/echickle/data/ZTF_Eclipses.txt'
    # out_dir   = sect_dir+'s%04d-lc-ZTF/'%sector
    # os.makedirs(out_dir, exist_ok=True)
    # run_lc_extraction(data_dir, out_dir, wd_cat, cam=cam, ccd=ccd,
    #                  mult_output=mult_output, tica=tica)

    wd_cat    = '/home/echickle/data/Gaia_blue.csv'
    out_dir   = sect_dir+'s%04d-lc-gaia/'%sector
    os.makedirs(out_dir, exist_ok=True)
    run_lc_extraction(data_dir, out_dir, wd_cat, cam=cam, ccd=ccd,
                     mult_output=mult_output, tica=tica)
    
    os.system('rm -r '+data_dir+'cam{}-ccd{}'.format(cam,ccd))

    
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
