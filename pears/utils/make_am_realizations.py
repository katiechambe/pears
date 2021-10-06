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

paths = SetupPaths()

# calculate the redshift at current snapshot
snapdata_path = f"{paths.path_snapdata}{sim}_snapdata.csv"
snapdata = QTable.read(snapdata_path)
redshift = snapdata['redshift'][snapdata['Snapshot'] == snapshot]

read_path = f"{paths.path_maxmass}{sim}_{phys}_{snapshot}.csv"

subhalo_data = QTable.read(read_path)



masses = np.array([1e10,1e12])

stels = AbundanceMatching(masses, redshift, 10).stellarMass()



# primID, primHaloMass, primStellars = [], [], []

haloGroup, haloID, haloMass, haloMaxMass, stellarMedian, stellarGauss = [], [], [], [], [], []

for ind in range(len(haloData)):
    haloGroup.append(haloData['Group Number'][ind])
    haloID.append(haloData['Subhalo ID'][ind])
    haloMass.append(haloData['Mass at snap'][ind])
    haloMaxMass.append(haloData['Max Mass'][ind]) # in 1e10 Msun

    maxMass         = haloData['Max Mass'][ind]*1e10 # halo mass in Msun (not 1e10!)
    maxMassRedshift = haloData['Redshift at Max Mass'][ind]
    median          = AM(maxMass, red, 'median').stellarMass()/1e10
    stellars        = AMGauss(maxMass, red, numRealizations).stellarMass()/1e10 # in 1e10Msun
    stellarMedian.append(median)
    stellarGauss.append(stellars)

    if ind%1000 == 0:
        print('on count ', ind, 'of ',len(haloData))



columnNames = ['Group Number','Subhalo ID', 'Mass at snap', 'Max Mass','Median Stellar'] + ['Stellar %i' % i for i in range(numRealizations)]
haloTable = np.column_stack([haloGroup, haloID, haloMass, haloMaxMass, stellarMedian, stellarGauss])
print('Created table')
haloDF = pd.DataFrame(data = haloTable, columns=columnNames)
print('Created dataframe')
print('Exporting to file')
haloDF.to_csv(fileToSave,index=False,header=True)
print('All done!')

save_path = f"{paths.path_am_mass}{sim}_{phys}_{snapshot}.csv"
