import numpy as np
import matplotlib.pyplot as plt
from astropy.io.votable import parse_single_table
import astropy.units as u
from astropy.constants import c
import pyphot
from pyphot import (unit, Filter)

# u,g,r,i,z
sky = np.array([22.96, 22.26, 21.2, 20.48, 19.6]) # from LSST
fwhm = np.array([0.75, 0.75, 0.75, 0.75, 0.75]) # estimate for ULTRACAM (=PSF=seeing)
# fwhm = np.array([0.75, 0.75, 0.75, 0.75, 0.75])*2 # estimate for ULTRACAM (=PSF=seeing)

zpt = np.array([24.15, 26.25, 25.7, 25.6, 24.77]) # vikdhillon NTT calculator

# pixel_size=0.35 on NTT, 0.3 on WHT
def SNR(source_mag,exp_time,n_frames,fwhm=fwhm,sky=sky,zpt=zpt,pixel_size=0.35,
        read_noise=3.5,dark_current=0.1):
    pixel_area=(pixel_size)**2
    sky_flux=2.512**(zpt-sky)*exp_time*n_frames
    source_flux=2.512**(zpt-source_mag)*exp_time*n_frames
    n_pixels=2.266*(fwhm/pixel_size)**2
    SNR=source_flux/np.sqrt(source_flux+sky_flux*pixel_area*n_pixels+\
                            read_noise**2*n_frames*n_pixels+\
                            dark_current*exp_time*n_frames)
    return SNR

# synthetic SDSS values from Gaia
votable = parse_single_table("/home/echickle/work/ATLASJ1138-5139/vizier_votable.vot")
sed_filter = votable.array['sed_filter']
sed_flux = votable.array['sed_flux']*u.Jy
sed_freq = votable.array['sed_freq']*u.gigahertz

SDSSu_freq = sed_freq[np.nonzero(sed_filter=='SDSS:u')][0]
SDSSg_freq = sed_freq[np.nonzero(sed_filter=='SDSS:g')][0]
SDSSr_freq = sed_freq[np.nonzero(sed_filter=='SDSS:r')][0]
SDSSi_freq = sed_freq[np.nonzero(sed_filter=='SDSS:i')][0]
SDSSz_freq = sed_freq[np.nonzero(sed_filter=='SDSS:z')][0]

SDSSu_wave = (c/SDSSu_freq).to(u.AA)
SDSSg_wave = (c/SDSSg_freq).to(u.AA)
SDSSr_wave = (c/SDSSr_freq).to(u.AA)
SDSSi_wave = (c/SDSSi_freq).to(u.AA)
SDSSz_wave = (c/SDSSz_freq).to(u.AA)

SDSSu_flux = sed_flux[np.nonzero(sed_filter=='SDSS:u')][0]
SDSSg_flux = sed_flux[np.nonzero(sed_filter=='SDSS:g')][0]
SDSSr_flux = sed_flux[np.nonzero(sed_filter=='SDSS:r')][0]
SDSSi_flux = sed_flux[np.nonzero(sed_filter=='SDSS:i')][0]
SDSSz_flux = sed_flux[np.nonzero(sed_filter=='SDSS:z')][0]

SDSSu_flux = SDSSu_flux.to(u.erg/u.cm**2/u.s/u.AA,
                           equivalencies=u.spectral_density(SDSSu_wave)).value
SDSSg_flux = SDSSg_flux.to(u.erg/u.cm**2/u.s/u.AA,
                           equivalencies=u.spectral_density(SDSSg_wave)).value
SDSSr_flux = SDSSr_flux.to(u.erg/u.cm**2/u.s/u.AA,
                           equivalencies=u.spectral_density(SDSSr_wave)).value
SDSSi_flux = SDSSi_flux.to(u.erg/u.cm**2/u.s/u.AA,
                           equivalencies=u.spectral_density(SDSSi_wave)).value
SDSSz_flux = SDSSz_flux.to(u.erg/u.cm**2/u.s/u.AA,
                           equivalencies=u.spectral_density(SDSSz_wave)).value

lib = pyphot.get_library()
SDSSu=lib['SDSS_u']
SDSSg=lib['SDSS_g']
SDSSr=lib['SDSS_r']
SDSSi=lib['SDSS_i']
SDSSz=lib['SDSS_z']

SDSSu_mag = -2.5 * np.log10(SDSSu_flux) - SDSSu.AB_zero_mag
SDSSg_mag = -2.5 * np.log10(SDSSg_flux) - SDSSg.AB_zero_mag
SDSSr_mag = -2.5 * np.log10(SDSSr_flux) - SDSSr.AB_zero_mag
SDSSi_mag = -2.5 * np.log10(SDSSi_flux) - SDSSi.AB_zero_mag
SDSSz_mag = -2.5 * np.log10(SDSSz_flux) - SDSSz.AB_zero_mag

mag=np.array([SDSSu_mag,SDSSg_mag,SDSSr_mag,SDSSi_mag,SDSSz_mag])
n_frames=1

# check counts match what we expect 2023_03_08
exp_time = 6.158917943838719
pixel_size=0.35
pixel_area=(pixel_size)**2
source_flux=2.512**(zpt-mag)*exp_time*n_frames
n_pixels=2.266*(fwhm/pixel_size)**2
print('Counts for 6.2s exposure: '+str(source_flux))

# print SNR for 3s exposure
exp_time=3
snr=SNR(mag,exp_time,n_frames)
print('SNR for 3s exposure: '+str(snr))

t=np.logspace(-1,1.25,100)

snr_read = []
snr_noread = []
for i in range(len(t)):
    exp_time = t[i]
    snr=SNR(mag,exp_time,n_frames)
    snr_read.append(snr)
    snr=SNR(mag,exp_time,n_frames,read_noise=0)
    snr_noread.append(snr)
snr_read = np.array(snr_read)
snr_noread = np.array(snr_noread)

plt.ion()
plt.figure()
plt.title('ULTRACAM SNR for source mag={}'.format(mag))
plt.plot(t, snr_noread[:,0], '-', label='No read noise')
plt.plot(t, snr_read[:,0], '-', label='Read noise=5e-')
plt.xlabel('Exposure time [seconds]')
plt.ylabel('SNR')
plt.yscale('log')
plt.xscale('log')
plt.legend()
plt.tight_layout()
plt.savefig('read_noise_floor_u.png', dpi=300)

plt.figure()
plt.title('ULTRACAM SNR for source mag={}'.format(mag))
plt.plot(t, snr_noread[:,1], '-', label='No read noise')
plt.plot(t, snr_read[:,1], '-', label='Read noise=5e-')
plt.yscale('log')
plt.xscale('log')
plt.xlabel('Exposure time [seconds]')
plt.ylabel('SNR')
plt.legend()
plt.tight_layout()
plt.savefig('read_noise_floor_g.png', dpi=300)

plt.figure()
plt.title('ULTRACAM SNR for source mag={}'.format(mag))
plt.plot(t, snr_noread[:,2], '-', label='No read noise')
plt.plot(t, snr_read[:,2], '-', label='Read noise=5e-')
plt.yscale('log')
plt.xscale('log')
plt.xlabel('Exposure time [seconds]')
plt.ylabel('SNR')
plt.legend()
plt.tight_layout()
plt.savefig('read_noise_floor_r.png', dpi=300)
