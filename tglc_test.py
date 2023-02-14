import os
from tglc.quick_lc import tglc_lc

target = "TIC 1201247611"
local_directory = "/data/submit/echickle/foo1/"
os.makedirs(local_directory, exist_ok=True)
tglc_lc(target=target, local_directory=local_directory, size=90, save_aper=False, limit_mag=20, get_all_lc=False, first_sector_only=False, sector=41, prior=None)
