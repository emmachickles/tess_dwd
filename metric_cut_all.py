import pdb
from vet_utils import *

data_dir = '/scratch/echickle/tess/BLS_results/'
sector_list = [56,57,58,59,60,61,62,63,64,65]

# Load Gaia white dwarf results 
result_list = append_result_file(data_dir, sector_list)

pdb.set_trace()
