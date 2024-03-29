import numpy as np
from astropy.coordinates import SkyCoord
import astropy.units as u
import pandas as pd
from astroquery.gaia import Gaia
import os

# data_dir = "/scratch/data/tess/lcur/ffi/"
# wd_tab = "/scratch/echickle/WDs.txt"
wd_tab = "/home/echickle/data/WDs.txt"

# >> Kevin's UCBs
ra1 = [234.883995, 86.6142665, 340.9290431, 66.45923447, 21.94843074, 84.5113413, 286.2972, 81.5434078, 83.38358801, 307.3429105, 110.5895062, 267.4803692, 337.1127, 296.5161979, 100.9031913, 100.0778794, 103.6232993, 322.7362856, 285.3559028, 357.8141279, 312.4636299, 273.2963257, 350.0851861, 313.8165708, 211.7342186, 259.5247867]

dec1 = [50.46077304, 38.7204002, 52.70166019, 38.98267612, 52.97027617, 19.8841512, 31.57565, 59.5792816, 2.153208645, 15.5752341, -18.6584643, 9.409030088, 49.82125, 32.05362787, 3.307603626, 17.6458415, 36.14832613, 44.34622882, 53.1581301, 63.0910333, 33.8647941, 42.86401367, 37.84185019, 46.8517723, 12.37868972, 34.1531371]

ra1, dec1 = np.array(ra1), np.array(dec1)
c= SkyCoord(ra=ra1*u.degree, dec=dec1*u.degree)

fnames = os.listdir('/data/ATLAS/')

# >> query all objects
Gaia.MAIN_GAIA_TABLE = "gaiadr3.gaia_source"
width = u.Quantity(0.1, u.deg)
height = u.Quantity(0.1, u.deg)
for i in range(len(c)):
    print(i)
    coord = c[i]
    r = Gaia.query_object_async(coordinate=coord, width=width, height=height)
    g = r['DESIGNATION'][0].split(' ')[-1]
    # r.pprint(max_lines=12, max_width=130)
    if g in fnames:
        print(c[i])
        print(g)


# wd_cat  = pd.read_csv(wd_tab, header=None, sep='\s+', dtype='str')
# ticid = wd_cat[0].to_numpy().astype('int')
# gid = wd_cat[3].to_numpy().astype('str')
# ra2, dec2 = wd_cat[1].to_numpy().astype('float'), wd_cat[2].to_numpy().astype('float')
# catalog= SkyCoord(ra=ra2*u.degree, dec=dec2*u.degree)

# idx, d2d, d3d = c.match_to_catalog_sky(catalog)
# matches= catalog[idx]
# dra, ddec = c.spherical_offsets_to(matches)
# sep = c.separation(matches)
# good_idx = np.nonzero(sep.to(u.arcsec).value < 2)
# print(matches[good_idx])
# print(ticid[idx][good_idx])
# print(gid[idx][good_idx])


# # mgab
# ra1 = [334.127233, 290.058818, 8.469317, 97.609023, 84.284387, 267.727448, 266.102959, 320.834172, 327.65423, 355.432715, 121.429093, 294.031619, 31.948969, 271.389582, 294.083066, 165.188166, 243.20429, 281.143375, 277.203223, 251.171596, 28.651189, 278.399813, 271.184634, 224.581366, 145.291371, 249.74499, 76.64014, 11.259771, 158.703887, 167.022366, 19.602117, 206.518421, 226.478111, 286.420908, 348.456363, 38.363497, 256.746138, 337.159227, 37.962464, 213.086238, 230.340651, 8.397302, 301.672242, 211.405559, 298.016194, 2.831444, 43.358164, 222.081, 125.046688, 248.587548, 217.472794, 150.06324, 126.265702, 307.048498, 320.331068, 101.684087, 29.67087, 62.569987, 121.174886, 118.229238, 319.081993, 93.14273, 156.722897, 135.096042, 210.574817, 61.220299, 15.510167, 114.61255, 261.349987, 194.312555, 150.399318, 5.861729, 98.125808, 323.051803, 132.489197, 264.486638, 298.793387, 200.620563, 7.694166, 61.649681, 145.2776, 289.005963, 96.643381, 297.303818, 78.390705, 336.128181, 101.011079, 357.121207, 15.309806, 105.13472, 293.209267, 0.654137, 226.907683]

# dec1 = [10.998309, 27.371746, 38.924915, 41.798482, -24.837474, 18.628736, 39.037842, 5.714922, 23.666681, 45.408752, -14.510147, 31.918566, 7.047877, 57.95545, 54.155878, 52.178829, 68.772309, 48.960087, 23.143996, 24.57459, 24.680022, 58.204948, 13.465076, 13.223986, 16.285595, 1.03998, 73.11303, 50.569015, 0.867048, 65.369847, 34.751205, 49.832578, -21.449036, 27.064795, 18.256609, 28.417863, 21.703238, 26.110494, 27.463674, 65.689804, 32.186628, 39.285167, 27.480695, 10.65535, 55.06785, -27.646677, 52.598155, -2.106623, 58.709781, -27.222424, 52.398122, 30.72515, 13.089111, 35.0903, 45.164227, 23.338327, -6.478306, -8.572031, -2.262536, -6.247477, 41.737626, -19.186594, -10.225072, 43.80374, -6.706418, 11.216442, -5.898116, 15.107347, 8.777635, 1.723506, -17.657502, -12.601798, 30.537165, 14.331333, 46.76276, 18.234318, 25.744812, 1.01466, -20.183986, -16.838307, 11.315382, -12.569772, 37.813012, 10.01346, -10.577674, 9.295982, 13.030307, -9.416927, 21.648836, 21.873519, 30.610314, 47.285249, 40.830757]

# ra1, dec1 = np.array(ra1), np.array(dec1)
# c= SkyCoord(ra=ra1*u.degree, dec=dec1*u.degree)

# for sector in [61, 62]:
#     for cam in [1,2,3,4]:
#         for ccd in [1,2,3,4]:
#             co_list = np.load(data_dir+"s%04d-lc/"%sector+"co-"+str(cam)+"-"+str(ccd)+".npy")
#             ticid   = np.load(data_dir+"s%04d-lc/"%sector+"id-"+str(cam)+"-"+str(ccd)+".npy")
#             ra2 = co_list[:,0]
#             dec2 = co_list[:,1]
#             catalog = SkyCoord(ra=ra2*u.degree, dec=dec2*u.degree)
#             idx, d2d, d3d = c.match_to_catalog_sky(catalog)
#             matches= catalog[idx]
#             dra, ddec = c.spherical_offsets_to(matches)
#             sep = c.separation(matches)
#             good_idx = np.nonzero(sep.to(u.arcsec).value < 2)
#             if len(good_idx[0]) > 0:
#                 print("S{}-{}-{}".format(sector, cam, ccd)+":")
#                 print(matches[good_idx])
#                 print(ticid[idx][good_idx])
