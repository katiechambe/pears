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


paths = SetupPaths()

read_path = f"{paths.path_maxmass}/{sim}_{phys}_{snapshot}.csv"
save_path = f"{paths.path_am_mass}/{sim}_{phys}_{snapshot}.csv"










snap = int(sys.argv[1]) # reading in the snapshot number
numRealizations = int(sys.argv[2])

#####################################
# defining the hardcoded paths, etc #
#####################################
baseDir      = '/xdisk/gbesla/mig2020/rsgrps/gbeslastudents/katie/dwarfPairsDark/'

spacingFile = baseDir + 'data/raw/data_spacing.csv'
primToRead = baseDir + 'data/processed/allPairs/primaries_'+str(snap)+'.csv'
secoToRead = baseDir + 'data/processed/allPairs/secondaries_'+str(snap)+'.csv'

primToSave = baseDir + 'data/processed/amRealizations/primaries_'+str(snap)+'.csv'
secoToSave = baseDir + 'data/processed/amRealizations/secondaries_'+str(snap)+'.csv'

#####################################################
# defining the snapshots and redshifts in Illustris #
#####################################################
snapData  = pd.read_csv(spacingFile)
snapshots = snapData['Snapshot'].values
scales    = snapData['Scale'].values
redshifts = snapData['Redshift'].values

primData = pd.read_csv(primToRead)
secoData = pd.read_csv(secoToRead)

primID, primHaloMass, primStellars = [], [], []
secoID, secoHaloMass, secoStellars = [], [], []

for ind in range(len(primData)):
    ## primaries first! ##
    primID.append(primData['Subhalo ID'][ind])
    primHaloMass.append(primData['Max Mass'][ind]) # in 1e10 Msun

    # to get the number of realizations I want, need array of masses to feed into AM class
    primMaxMass         = primData['Max Mass'][ind]*1e10 # halo mass in Msun (not 1e10!)
    primMaxMassRedshift = primData['Redshift at Max Mass'][ind]
    primStellarMasses   = np.array([(AM(primMaxMass, primMaxMassRedshift, 'gauss').stellarMass())[0]/1e10 for i in range(numRealizations)])
    primStellars.append(primStellarMasses)

    ## secondaries second! ##
    secoID.append(secoData['Subhalo ID'][ind])
    secoHaloMass.append(secoData['Max Mass'][ind]) # in 1e10 Msun

    # to get the number of realizations I want, need array of masses to feed into AM class
    secoMaxMass         = secoData['Max Mass'][ind]*1e10 # halo mass in Msun (not 1e10!)
    secoMaxMassRedshift = secoData['Redshift at Max Mass'][ind]
    secoStellarMasses   = np.array([(AM(secoMaxMass, secoMaxMassRedshift, 'gauss').stellarMass())[0]/1e10 for i in range(numRealizations)])
    secoStellars.append(secoStellarMasses)

columnNames = ['Subhalo ID', 'Max Mass'] + ['Stellar %i' % i for i in range(numRealizations)]

primTable = np.column_stack([primID, primHaloMass, primStellars])
primDF = pd.DataFrame(data = primTable, columns=columnNames)
primDF.to_csv(primToSave,index=False,header=True)

secoTable = np.column_stack([secoID, secoHaloMass, secoStellars])
secoDF = pd.DataFrame(data = secoTable, columns=columnNames)
secoDF.to_csv(secoToSave,index=False,header=True)
