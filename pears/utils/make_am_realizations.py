""" 
Creates stellar mass realizations for DM halos

Dependencies:
-------------
Requires max mass files to exist! 
Use make_max_mass_file.py before running this script

"""
__author__ = "Katie Chamberlain"
__status__ = "Beta - forever~"
__date__   = "October 2021"

import sys
import h5py
from astropy.table import QTable
from utils.abundance_matching import AbundanceMatching
from utils.paths import SetupPaths

snapshot = int(sys.argv[1])
sim = str(sys.argv[2])
physics = str(sys.argv[3])
num_reals = int(sys.argv[4])

paths = SetupPaths()

# importing paths
snapdata_path = f"{paths.path_snapdata}{sim}_snapdata.csv"
read_path = f"{paths.path_maxmass}{sim}_{physics}_{snapshot}.csv"

# defining redshift corresponding to snapshot
snapdata = QTable.read(snapdata_path)
redshift = snapdata['redshift'][snapdata['snapshot'] == snapshot][0]

# reading in subhalo information
subhalo_data = QTable.read(read_path)

ids, masses, maxmasses, maxsnaps, stellar_reals = [], [], [], [], []

for ind in range(len(subhalo_data)):
    halo_id = subhalo_data['Subhalo ID'][ind]
    halo_mass = subhalo_data['Current Snap Mass'][ind]
    halo_maxmass = subhalo_data['Max Mass'][ind]
    halo_snap = subhalo_data['Max Mass Snap'][ind]

    ids.append(halo_id)
    masses.append(halo_mass)
    maxmasses.append(halo_maxmass)
    maxsnaps.append(halo_snap)
        
    stars = AbundanceMatching(maxmass=halo_maxmass*1e10, 
                              redshift=redshift,
                              samples=num_reals).stellar_mass()[0] # in 1e10Msun

    stellar_reals.append(stars/1e10)

    if ind%10000 == 0:
        print('on count ', ind, 'of ',len(subhalo_data))

data_dict = {"Subhalo ID":ids, 
             "Max Mass":maxmasses, 
             "Max Mass Snap":maxsnaps, 
             "Current Snap Mass":masses, 
             "Stellar Masses":stellar_reals}

save_path = f"{paths.path_am_mass}{sim}_{physics}_{snapshot}.hdf5"

f = h5py.File(save_path, 'w')
for key, val in data_dict.items():
    f[key] = val
f.close()

print(f"Saved at {paths.path_am_mass}{sim}_{physics}_{snapshot}.hdf5")
