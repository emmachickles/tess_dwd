import numpy as np
import matplotlib.pyplot as plt
import ellc

period = 26.  # Orbital period in minutes
mass_ratio = 0.3 # Mass ratio
inclination = 80. # Inclination angle in degrees
radius_1 = 0.1 # Radius of star 1 (in units of separation) 
radius_2 = 0.2 # Radius of star 2 (in units of separation)
sbratio = 0.6 # Surface brightness ratio 

fig, ax = plt.subplots()

# time=np.random.uniform(0,80,1000)
time=np.linspace(0,80,10000)

# Set up the ELLC model
flux = ellc.lc(
    t_obs=time,  # Observed time array                
    period=period,  # Orbital period in minutes
    radius_1=radius_1,  # Radius of star 1 (in units of separation)
    radius_2=radius_2,  # Radius of star 2 (in units of separation)
    sbratio=sbratio,  # Surface brightness ratio
    incl=inclination,  # Inclination angle in degrees
    q=mass_ratio,  # Mass ratio
    f_c=0.2,  # Limb darkening coefficients for star 1 (white dwarf)
    f_s=0.6,  # Limb darkening coefficients for star 2 (secondary)
    gdc_1=0.1,  # Gravity darkening coefficients for star 1 (white dwarf)
    gdc_2=0.8  # Gravity darkening coefficients for star 2 (secondary)
)

# Plot the light curve
ax.plot(time, flux, '.', ms=1)
ax.set_title(f"Period: {period} min, Mass ratio: {mass_ratio}")
ax.set_xlabel('Time [minutes]')
ax.set_ylabel('Flux')
fig.savefig('/home/echickle/foo.png')
