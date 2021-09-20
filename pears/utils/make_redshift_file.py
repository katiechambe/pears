""" 
Collects the redshift, snapshot, and scale data for the simulation

Dependencies:
-------------
ReadCats class from read_group_cats.py
"""
__author__ = "Katie Chamberlain"
__status__ = "Beta - forever~"
__date__   = "September 2021"

import numpy as np
from utils.read_group_cats import ReadCats
from astropy.table import QTable

snapshot_nums = np.arange(0,135)

reds = []
snapshots = []
scales = [] 

sim = "Illustris"

for snap in snapshot_nums:
    try:
        redshift = ReadCats(snapshot=snap).redshift
        scale = 1/(1+redshift)

        reds.append(redshift)
        scales.append(scale)
        snapshots.append(snap)
    except :
        continue

t = QTable()
t['snapshot'] = snapshots
t['redshift'] = reds
t['scale'] = scales
t.write("../../data/.csv"
        overwrite=True)

print(f"Saved data")
