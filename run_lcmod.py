import ellc
import numpy as np
import matplotlib.pyplot as plt

def plot_lc(t, y, output_dir, title='', suffix=''):
    plt.figure()
    plt.plot(t, y, '.k', ms=1)
    plt.title(title)
    plt.xlabel('Time [s]')
    plt.ylabel('Relative flux')
    fname = output_dir+'lc'+suffix+'.png'
    plt.savefig(fname, dpi=300)
    print('Saved '+fname)


output_dir = '/scratch/echickle/dwd/'

# based off https://github.com/pmaxted/ellc/blob/master/examples/batman/run_batman.py

# a       = 10.                 # >> semi-major axis (in units of stellar radii)
# r_1     = 1/a
# rp      = 0.1                 # >> planet radious (in units of stellar radius)
# r_2     = r_1*rp
# incl    = 90.              
# ecc     = 0.1                 # >> eccentricity
# w       = 69.                 # >> longitude of periastron (in degrees) 
# sbratio = 0                   # >> surface brightness ratio, S_2/S_1
# ld_1    = 'quad'              # >> limb darkening model
# ldc_1   = [0.1, 0.3]          # >> limb darkening coefficients

# t   = np.linspace(-0.05, 0.05, 1000)
# f_s = np.sqrt(ecc) * np.sin(w * np.pi / 180.)
# f_c = np.sqrt(ecc) * np.cos(w * np.pi / 180.) 

# flux = ellc.lc(t,radius_1=r_1, radius_2=r_2,incl=incl,sbratio=sbratio,
#                ld_1=ld_1, ldc_1=ldc_1,shape_1='sphere',shape_2='sphere',
#                grid_1='default',grid_2='default', f_s=f_s, f_c=f_c)

# plt.figure()
# plt.plot(t, flux, '.k')
# plt.title('r2=0.1r1')
# plt.xlabel('Time')
# plt.ylabel('Flux')
# plt.savefig(output_dir+'lc1.png')

# r_2     = r_1

# flux = ellc.lc(t,radius_1=r_1, radius_2=r_2,incl=incl,sbratio=sbratio,
#                ld_1=ld_1, ldc_1=ldc_1,shape_1='sphere',shape_2='sphere',
#                grid_1='default',grid_2='default', f_s=f_s, f_c=f_c)

# plt.figure()
# plt.plot(t, flux, '.k')
# plt.title('r2=1r2')
# plt.xlabel('Time')
# plt.ylabel('Flux')
# plt.savefig(output_dir+'lc2.png')


# r_2     = r_1 * 0.7

# flux = ellc.lc(t,radius_1=r_1, radius_2=r_2,incl=incl,sbratio=sbratio,
#                ld_1=ld_1, ldc_1=ldc_1,shape_1='sphere',shape_2='sphere',
#                grid_1='default',grid_2='default', f_s=f_s, f_c=f_c)

# plt.figure()
# plt.plot(t, flux, '.k')
# plt.title('r2=0.7r2')
# plt.xlabel('Time')
# plt.ylabel('Flux')
# plt.savefig(output_dir+'lc3.png')

# a = 11.218e-2 # >> in solar radii
# r_1 = 1.562e-2 / a # >> in units of a
# r_2 = 3.140e-2 / a
# sbratio = (10000 / 48900) ** 4
# incl = 84.15 # >> deg
# period = 414.79 # >> s
# t   = np.linspace(-500, 500, 10000) # >> s
# m_1 = 0.610 # >> in solar masses
# m_2 = 0.210
# q = m_2/m_1
# heat_2 = 1 # >> =A_g/2 where A_g is the geometric albedo


# flux = ellc.lc(t,radius_1=r_1, radius_2=r_2,incl=incl,sbratio=sbratio,
#                ld_1=ld_1, ldc_1=ldc_1,shape_1='roche',shape_2='roche',
#                grid_1='default',grid_2='default', f_s=f_s, f_c=f_c,
#                period=period, a=a, q=q, heat_2=heat_2)
# plot_lc(t, flux, output_dir, title='7 min', suffix='_7min')

a = 0.1227 # >> in solar radii
r_1 = 0.0298 / a # >> in units of a
r_2 = 0.0275 / a
sbratio = (19200 / 26300)**4 # >> s_2 / s_1
incl = 82.12 # >> deg
period = 527.93 # >> s
# t   = np.linspace(-1000, 1000,20000) # >> s
t = np.linspace(-period, period, 1000)
m_1 = 0.323 # >> in solar masses
m_2 = 0.335
q = m_2/m_1
heat_2 = 0.5


flux = ellc.lc(t,radius_1=r_1, radius_2=r_2,incl=incl,sbratio=sbratio,
               ld_1=ld_1, ldc_1=ldc_1,shape_1='roche',shape_2='roche',
               grid_1='default',grid_2='default', f_s=f_s, f_c=f_c,
               period=period, a=a, q=q, heat_2=heat_2)
plot_lc(t, flux, output_dir, title='8.8 min', suffix='_8min')
