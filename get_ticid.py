import pandas as pd
from astropy.coordinates import SkyCoord
from astropy import units as u
from astroquery.mast import Catalogs

df = pd.read_csv('/home/echickle/data/Gaia_blue.csv')

def query_ticids_for_coordinates(ra, dec):
    result = Catalogs.query_region('{} {}'.format(ra, dec), catalog="Tic", radius=1 * u.arcsec)
    tic_ids = result["ID"].tolist()

    return tic_ids

# Example usage
ra = 123.456
dec = 12.345
tic_ids = query_ticids_for_coordinates(ra, dec)
print(tic_ids)
