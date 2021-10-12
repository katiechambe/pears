""" 
Creates stellar mass realizations for DM halos

Dependencies:
-------------
ReadCats class from read_group_cats.py

"""
__author__ = "Katie Chamberlain"
__status__ = "Beta - forever~"
__date__   = "October 2021"

import numpy as np
import pandas as pd
import sys
import astropy.units as u
from astropy.table import QTable
from utils.abundance_matching import AbundanceMatching
from utils.paths import SetupPaths


#TODO: define physics, simulation, and snapshot
#TODO: allow for median in abundance matching code! as alt calc.
sim = "Illustris" 
physics = "dark"
snapshot = 135
num_reals = 1000

paths = SetupPaths()

# calculate the redshift at current snapshot
snapdata_path = f"{paths.path_snapdata}{sim}_snapdata.csv"
snapdata = QTable.read(snapdata_path)
redshift = snapdata['redshift'][snapdata['snapshot'] == snapshot][0]

read_path = f"{paths.path_maxmass}{sim}_{physics}_{snapshot}.csv"

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
        
        # median          = AM(maxMass, red, 'median').stellarMass()/1e10
    stars = AbundanceMatching(maxmass=halo_maxmass*1e10, 
                              redshift=redshift,
                              samples=num_reals).stellar_mass()[0] # in 1e10Msun

    # NOTE: am i going to have to reshape this? what type does this come out as? 
    # What are the dimensions of the data structure that comes out from this.
    # Is it easier to save the stellar masses as one big stucture and then reshape it? 
    stellar_reals.append(stars/1e10)

    if ind%5000 == 0:
        print('on count ', ind, 'of ',len(subhalo_data))

t = QTable()
t["Subhalo ID"] = ids
t["Max Mass"] = maxmasses * u.Unit(1e10 * u.Msun)
t["Max Mass Snap"] = maxsnaps
t["Current Snap Mass"] = masses * u.Unit(1e10 * u.Msun)
t["Stellar Masses"] = stellar_reals * u.Unit(1e10 * u.Msun)
save_path = f"{paths.path_am_mass}{sim}_{physics}_{snapshot}.ecsv"
t.write(save_path, overwrite=True)
print('All done!')

