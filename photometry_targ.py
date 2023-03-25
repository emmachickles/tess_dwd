import os

from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy import wcs
from astropy.io import fits

from photutils import aperture_photometry, CircularAperture, SkyCircularAnnulus
from photutils import SkyCircularAperture

ra = 228.3879351 # >> can also be a list of radecs in the same camera
dec = 70.6227716
out_dir = '/data/submit/echickle/data/' # >> where lc, ffi directory will be saved

# >> look at https://tess.mit.edu/observations/
sector = 41
cam = 3

def download_cam_ffi(out_dir, sector, cam):
    '''Note that each FFI file is ~0.034GB, and each camera > 100 GB.'''
    
    # >> get cURL File for Full Frame Images
    curl_f = 'tesscurl_sector_'+str(sector)+'_ffic.sh'
    os.system('curl -O https://archive.stsci.edu/missions/tess/'+\
              'download_scripts/sector/'+curl_f)

    # >> download FFIs in camera
    with open(curl_f, 'r') as f:
        lines = f.readlines()[1:]
    cam_dir = 's'+str(sector)+'-cam'+str(cam)+'/'
    if not os.path.exists(out_dir+cam_dir):
        os.makedirs(out_dir+cam_dir)
    for line in lines:
        ffi_cam = int(line.split(' ')[5].split('-')[2])
        if ffi_cam == cam: 
            line = line.split(' ')
            line[5] = out_dir+cam_dir+line[5]
            line = ' '.join(line)
            os.system(line)

def extract_lc(out_dir, sector, cam, ra, dec):
    # >> star coordinate
    coord = SkyCoord(ra=ra*u.degree, dec=dec*u.degree, frame='icrs')

    # >> radii of aperture and background annulus
    N_ap=0.7
    aperture = SkyCircularAperture(trimmed_catalog, r=N_ap*21 * u.arcsec)
    N_in=1.5
    N_out=2
    background = SkyCircularAnnulus(trimmed_catalog,
                                    N_in*21*u.arcsec,N_out*21*u.arcsec)

    # >> get ffi file names
    cam_dir = 's'+str(sector)+'-cam'+str(cam)+'/'
    ffi_fnames = os.listdir(out_dir+cam_dir)
    
    # >> get flux
    time, flux = [], []
    for f in ffi_fnames:
        try:
            hdul = fits.open(f)            
            t = hdul[0].header['TSTART']+300/86400.0 # observation time in BTJD
            w = wcs.WCS(hdu_list[1].header)

            image = hdu_list[1].data
            n_dty = dt.shape[0]
            n_dtx = dt.shape[1]

            pix_aperture = aperture.to_pixel(w)
            phot_table = aperture_photometry(image, pix_aperture)  

            background_pix_aperture = background.to_pixel(w)
            background_phot_table = aperture_photometry(image, background_pix_aperture)  

            norm=background_pix_aperture.area/pix_aperture.area

            for col in phot_table.colnames:  
                phot_table[col].info.format = '%.8g'  # for consistent table out

            y = phot_table['aperture_sum'].value-background_phot_table['aperture_sum'].value/norm
            
            time.append(t)
            flux.append(y)


            import pdb
            pdb.set_trace()
        except:
            pass

    # -- save timeseries -------------------------------------------------------
    time = np.array(time)
    flux = np.array(flux)
    np.save(out_dir+'lc_'+str(ra)+'_'+str(dec)+'.npy', np.array([time, flux]).T)
    return

# download_cam_ffi(out_dir, sector, cam)
extract_lc(out_dir, sector, cam, ra, dec)
