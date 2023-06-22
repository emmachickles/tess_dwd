import numpy as np
import matplotlib.pyplot as plt
import ellc
import pdb

# Define the range of parameters to explore
periods = [12., 26., 61.]  # Orbital periods in minutes
mass_ratios = [0.3, 0.5]  # Mass ratios (secondary mass / primary mass)
inclinations = [90., 85., 80.]  # Inclination angles in degrees

q = np.array(mass_ratios)
rl = 0.49*q**(2/3) / (0.6*q**(2/3) + np.log(1+q**(1/3))) # Roche lobe radius of secondary in units of separation
print(rl)

# Generate synthetic light curves for each combination of parameters
fig, axes = plt.subplots(len(periods), len(mass_ratios), figsize=(12, 9), sharex=True, sharey=True)

for i, period in enumerate(periods):
    for j, mass_ratio in enumerate(mass_ratios):
        for k, inclination in enumerate(inclinations):
            # time=np.random.uniform(0,80,1000)
            time=np.linspace(0,80,10000)

            # Set up the ELLC model
            flux = ellc.lc(
                t_obs=time,  # Observed time array                
                period=period,  # Orbital period in minutes
                radius_1=0.1,  # Radius of star 1 (in units of separation)
                radius_2=0.2,  # Radius of star 2 (in units of separation)
                sbratio=0.6,  # Surface brightness ratio
                incl=inclination,  # Inclination angle in degrees
                q=mass_ratio,  # Mass ratio
                f_c=0.2,  # Limb darkening coefficients for star 1 (white dwarf)
                f_s=0.6,  # Limb darkening coefficients for star 2 (secondary)
                gdc_1=0.1,  # Gravity darkening coefficients for star 1 (white dwarf)
                gdc_2=0.8  # Gravity darkening coefficients for star 2 (secondary)
            )

            # Plot the light curve
            ax = axes[i, j]
            ax.plot(time, flux, '.', ms=1, label=f"Inc: {inclination}°")
            ax.set_title(f"Period: {period} min, Mass ratio: {mass_ratio}")

# Set common labels and adjust plot spacing
fig.text(0.4, 0.04, 'Time', ha='center')
fig.text(0.04, 0.5, 'Flux', va='center', rotation='vertical')
plt.subplots_adjust(wspace=0.1, hspace=0.3)

# Add legend and display the plot
plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
plt.savefig('/home/echickle/foo.png')
pdb.set_trace()
