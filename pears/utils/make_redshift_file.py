""" 
Collects the redshift, snapshot, and scale data for the simulation

Dependencies:
-------------
ReadCats class from read_group_cats.py
"""
__author__ = "Katie Chamberlain"
__status__ = "Beta - forever~"
__date__   = "September 2021"

import utils.readsubfHDF5Py3 as readSub
from astropy.table import QTable
import numpy as np
from utils.paths import SetupPaths



reds = []
snapshots = []
scales = [] 

paths = SetupPaths()

# sim = "Illustris"
sim = "TNG"

if sim == "Illustris":
    catpath = paths.path_illustrisdark
    backup = paths.path_illustrishydro
    snapshot_nums = np.arange(0,136)

        
elif sim == "TNG":
    catpath = paths.path_tngdark
    backup = paths.path_tnghydro
    snapshot_nums = np.arange(0,100)


for snap in snapshot_nums:
    print(snap)
    try:
        catalog = readSub.subfind_catalog(
        basedir=catpath,
        snapnum=snap,
        keysel=[]
        )

        redshift = catalog.redshift
        scale = 1/(1+redshift)

        reds.append(redshift)
        scales.append(scale)
        snapshots.append(snap)
    except OSError:
        print("snapshot corrupt in dark")
        try:
            catalog = readSub.subfind_catalog(
            basedir=backup,
            snapnum=snap,
            keysel=[]
            )

            redshift = catalog.redshift
            scale = 1/(1+redshift)

            reds.append(redshift)
            scales.append(scale)
            snapshots.append(snap)
        except OSError:
            print(f"snapshot {snap} does not exist")
            continue
        except SystemExit:
            continue

t = QTable()
t['snapshot'] = snapshots
t['redshift'] = reds
t['scale'] = scales
t.write(f"{paths.path_snapdata}{sim}_snapdata.csv",
        overwrite=True)

print(f"Saved data")


