import os
from tglc.quick_lc import tglc_lc

local_directory = "/data/submit/echickle/kbucb/"
target = "TIC 2040677137"
tglc_lc(target=target, local_directory=local_directory, size=90,
        save_aper=False, limit_mag=20, get_all_lc=False,
        first_sector_only=False, sector=57,
        prior=None)
