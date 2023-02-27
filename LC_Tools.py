import numpy as np
from astropy.time import Time
from astropy.coordinates import EarthLocation
from astropy.coordinates import SkyCoord  # High-level coordinates
from astropy.coordinates import ICRS, Galactic, FK4, FK5, BarycentricTrueEcliptic  # Low-level frames
from astropy.coordinates import Angle, Latitude, Longitude  # Angles
import astropy.units as u


"Function for converting a series of timestamps to Barycentric Julian Date format in Barycentric Dynamical time"

def BJDConvert(times, RA, Dec, date_format='mjd', telescope='Palomar'):

	t = Time(times,format=date_format,scale='utc')

	t2=t.tcb

	c = SkyCoord(RA,Dec, unit="deg")

	d=c.transform_to(BarycentricTrueEcliptic)

	Observatory=EarthLocation.of_site(telescope)

	delta=t2.light_travel_time(c,kind='barycentric',location=Observatory)

	BJD_TCB=t2+delta

	return BJD_TCB

"Function which returns phases corresponding to timestamps in a lightcurve given a period P, period derivative Pdot, and reference epoch t0"

def PhaseFold(times, P, Pdot=0, t0=0):


	"If no reference epoch is supplied, reference epoch is set to earliest time in lightcurve"
	if t0==0:
		times=times-np.min(times)
	else:
		times=times-t0
	phases=((times-1/2*Pdot/P*(times)**2) % P)/P

	return phases

"Weighted binning of a lightcurve with times t, fluxes y, and flux errors dy, at period p with N bins."

def binning(t,y,dy,P,Pdot=0,t0=0,N=500,cycles=3,Norm=False):
    binned_LC=[]
    mean_phases=np.linspace(0,1-1/N,N)
    phases=PhaseFold(t, P, Pdot, t0)
    lightcurve=np.array((phases,y,dy)).T
    lightcurve=lightcurve[np.argsort(lightcurve[:,0])]
    for i in mean_phases:
        lightcurve_bin=lightcurve[lightcurve[:,0]>i]
        lightcurve_bin=lightcurve_bin[lightcurve_bin[:,0]<i+1/N]
        weights=1/(lightcurve_bin[:,2]**2)
        weighted_mean_flux=np.sum(lightcurve_bin[:,1]*weights)/np.sum(weights)
        weighted_mean_flux_error=np.sqrt(1/np.sum(weights))
        binned_LC.append((i+0.5/N,weighted_mean_flux,weighted_mean_flux_error))
    binned_LC=np.array(binned_LC)
    if Norm=='Median':
    	binned_LC[:,2]=binned_LC[:,2]/np.median(binned_LC[:,1])
    	binned_LC[:,1]=binned_LC[:,1]/np.median(binned_LC[:,1])
    elif Norm=='Min':
    	binned_LC[:,2]=binned_LC[:,2]/np.min(binned_LC[:,1])
    	binned_LC[:,1]=binned_LC[:,1]/np.min(binned_LC[:,1])
    elif Norm=='Max':
    	binned_LC[:,2]=binned_LC[:,2]/np.max(binned_LC[:,1])
    	binned_LC[:,1]=binned_LC[:,1]/np.max(binned_LC[:,1])
    elif Norm=='Mean':
    	binned_LC[:,2]=binned_LC[:,2]/np.mean(binned_LC[:,1])
    	binned_LC[:,1]=binned_LC[:,1]/np.mean(binned_LC[:,1])
    	
    if cycles==1:
    	binned_LC=binned_LC
    elif cycles==2:
    	binned_LC2=np.array((binned_LC[:,0]-1,binned_LC[:,1],binned_LC[:,2])).T
    	binned_LC=np.hstack((binned_LC2,binned_LC))
    elif cycles==3:
    	binned_LC2=np.array((binned_LC[:,0]-1,binned_LC[:,1],binned_LC[:,2])).T
    	binned_LC3=np.array((binned_LC[:,0]+1,binned_LC[:,1],binned_LC[:,2])).T
    	binned_LC=np.vstack((binned_LC2,binned_LC))    
    	binned_LC=np.vstack((binned_LC,binned_LC3))  	

    return binned_LC
