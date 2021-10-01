""" 
Creates data file of maximum mass of all subhalos with mass >1e9 and 


--------
How to run:
-----------
python make_max_mass_file.py snapshot# simulation

where snapshot can range between 0 and 136 for "Illustris"
and 

"""
__author__ = "Katie Chamberlain"
__status__ = "Beta - forever~"
__date__   = "September 2021"

import numpy as np
import h5py
import sys
import os
from utils.merger_trees import TraceMergerTree
from utils.get_groups import GetGroups

snap = int(sys.argv[1])
sim = str(sys.argv[2])
# snap = 135
# sim = "Illustris"

if sim == "Illustris":
    little_h = 0.704

elif sim == "TNG":
    little_h = 0.6774

kwargs = {"group_mass_min": 8,
          "group_mass_max": 500,
          "little_h": little_h}

for phys in ["dark","hydro"]:
    groups = GetGroups(snapshot=snap,
                        sim=sim,
                        physics=phys,
                        **kwargs)

    subids = groups.subhalo_ids
    max_mass = []
    current_mass = []
    max_mass_snapshot = []


    for subid in subids[0:10]:
        max_mass_info = TraceMergerTree(snapshot=snap,
                            subfindID=subid,
                            sim=sim, 
                            physics=phys,
                            **kwargs)
        max_mass.append(max_mass_info.maxmass[0])
        max_mass_snapshot.append(max_mass_info.maxmass[1])

        current_mass.append(max_mass_info.masses[-1])

    file_name = f"{groups.path_data}max_masses/{sim}_{phys}_{snap}.hdf5"
    if os.path.exists(file_name):
        os.remove(file_name)

    savefile = h5py.File(file_name, "w")
    subs = savefile.create_dataset("Subhalo ID", 
                                   dtype=np.int64, 
                                   data=subids)
    maxs = savefile.create_dataset("Max Mass", 
                                   dtype=np.float64, 
                                   data=max_mass)
    snaps = savefile.create_dataset("Max Mass Snap", 
                                    dtype=np.int64, 
                                    data=max_mass_snapshot)
    masss = savefile.create_dataset("Current Snap Mass", 
                                    dtype=np.float64, 
                                    data=current_mass)

    for i in [maxs,masss]:
        i.attrs['units'] = "1e10 Msun"
    savefile.close()

    