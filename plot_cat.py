import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import astropy.coordinates as coord
import astropy.units as u

data_dir = '/data/submit/echickle/'
out_dir = '/home/submit/echickle/'

fig = plt.figure(figsize=(8,6))
ax = fig.add_subplot(111, projection="mollweide")
ax.grid(True)

# -- Gaia EDR3 WD Candidates ---------------------------------------------------
wd_cat  = pd.read_csv(data_dir+'WDs.txt', header=None, sep='\s+')
ra1 = wd_cat[1].to_numpy()
inds = np.nonzero(ra1 > 180.)
ra1[inds] = ra1[inds]-360

ra1 = coord.Angle(ra1*u.degree).radian
dec1 = wd_cat[2].to_numpy()
dec1 = coord.Angle(dec1*u.degree).radian

ax.plot(ra1, dec1, '.', ms=0.01, c='gray', label='Gaia eDR3 WD candidates')

# -- SDSS WD main-sequence binaries : J/MNRAS/382/1377 -------------------------
# JHHMMSS.ss+DDMMSS.s
# 'SDSSJ234459.62+002749.9'
with open(data_dir+'J_MNRAS_382_1377/table2.dat', 'r') as f:
    lines = f.readlines()
with open(data_dir+'J_MNRAS_458_3808/tables1.dat', 'r') as f:
    lines.extend(f.readlines())
ra2 = []
dec2 = []
for i in range(len(lines)):
    ra = lines[i].split(' ')[0][5:14]
    ra = ra[:2]+'h'+ra[2:4]+'m'+ra[4:]+'s'
    dec = lines[i].split(' ')[0][14:]
    dec = dec[:3]+'d'+dec[3:5]+'m'+dec[5:]+'s'
    co = coord.SkyCoord(ra, dec)
    # co = coord.SkyCoord.from_name(lines[i].split(' ')[0])
    if co.ra.value > 180.:
        ra2.append( (co.ra - 360*u.deg).radian)
    else:
        ra2.append(co.ra.radian)
    dec2.append(co.dec.radian)
ra2, dec2 = np.array(ra2), np.array(dec2)
ax.plot(ra2, dec2, '.', ms=1, label='SDSS WDMS binaries')

# -- SDSS AM CVn ---------------------------------------------------------------
with open(data_dir+'J_MNRAS_429_2143/tablea1.dat', 'r') as f:
    lines = f.readlines()
ra3 = []
dec3 = []
for i in range(len(lines)):
    ra = lines[i][:11]
    dec = lines[i][12:23]
    co = coord.SkyCoord(ra, dec, unit=(u.hourangle, u.deg))
    if co.ra.value > 180.:
        ra3.append( (co.ra - 360*u.deg).radian)
    else:
        ra3.append(co.ra.radian)
    dec3.append(co.dec.radian)
ra3, dec3 = np.array(ra3), np.array(dec3)
ax.plot(ra3, dec3, '.', ms=1, label='SDSS AM CVns')

# ------------------------------------------------------------------------------

for j in range(2):
    if j == 0: # >> sector 56
        s=56
        cam_ra = [345.364, 333.052, 311.763, 266.081]
        cam_dec = [14.958, 36.296, 55.430, 65.859]
        c='r'
    if j == 1:
        s = 57
        cam_ra = [11.279, 356.471, 325.102, 267.376]
        cam_dec = [26.108, 47.064, 63.839, 65.205]
        c = 'orange'

    for i in range(4):
        ra_rng = np.linspace(cam_ra[i] - 12, cam_ra[i] + 12, 100)
        inds = np.nonzero(ra_rng > 180.)
        ra_rng[inds] = ra_rng[inds]-360
        ra_rng = coord.Angle(ra_rng*u.degree).radian

        dec_rng = np.linspace(cam_dec[i] - 12, cam_dec[i] + 12, 100)
        dec_rng = coord.Angle(dec_rng*u.degree).radian

        ax.plot(ra_rng, np.ones(100)*np.min(dec_rng), '-', c=c)
        ax.plot(ra_rng, np.ones(100)*np.max(dec_rng), '-', c=c)

        ax.plot(np.ones(100)*np.min(ra_rng), dec_rng, '-', c=c)
        ax.plot(np.ones(100)*np.max(ra_rng), dec_rng, '-', c=c)    

    ax.legend(loc='upper right')
    plt.savefig(out_dir+'map.png', dpi=100)

