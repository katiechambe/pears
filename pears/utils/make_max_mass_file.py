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
from astropy.table import QTable
import astropy.units as u

snap = int(sys.argv[1])
sim = str(sys.argv[2])

if sim == "Illustris":
    little_h = 0.704

elif sim == "TNG":
    little_h = 0.6774

kwargs = {"group_mass_min": 8,
          "group_mass_max": 500,
          "little_h": little_h}

for phys in ["hydro","dark"]:
    groups = GetGroups(snapshot=snap,
                        sim=sim,
                        physics=phys,
                        **kwargs)

    subids = groups.subhalo_ids
    max_mass = []
    current_mass = []
    max_mass_snapshot = []

    print(len(subids))
    countertho = 0
    for subid in subids:
        countertho += 1
        if countertho%10000 == 0:
            print(f"Done with {countertho} out of{len(subids)}")
        try:
            # max mass and corresponding snapnumber
            max_mass_info = TraceMergerTree(snapshot=snap,
                            subfindID=subid,
                            sim=sim, 
                            physics=phys,
                            **kwargs)
    
        except AttributeError:
            print('Could not find merger tree for subhalo:', subid)
            print('Appending many zeros')
            max_mass.append(0)
            max_mass_snapshot.append(0)
            current_mass.append(0)

        else:
            max_mass.append(max_mass_info.maxmass[0])
            max_mass_snapshot.append(max_mass_info.maxmass[1])
            current_mass.append(max_mass_info.masses[0])

    file_name = f"{groups.path_data}max_masses/{sim}_{phys}_{snap}.ecsv"
    if os.path.exists(file_name):
        os.remove(file_name)

    # saving file as hdf5
    # savefile = h5py.File(file_name, "w")
    # subs = savefile.create_dataset("Subhalo ID", 
    #                                dtype=np.int64, 
    #                                data=subids)
    # maxs = savefile.create_dataset("Max Mass", 
    #                                dtype=np.float64, 
    #                                data=max_mass)
    # snaps = savefile.create_dataset("Max Mass Snap", 
    #                                 dtype=np.int64, 
    #                                 data=max_mass_snapshot)
    # masss = savefile.create_dataset("Current Snap Mass", 
    #                                 dtype=np.float64, 
    #                                 data=current_mass)

    # for i in [maxs,masss]:
    #     i.attrs['units'] = "1e10 Msun"
    # savefile.close()

    # saving file as ecsv
    t = QTable()
    t["Subhalo ID"] = subids
    t["Max Mass"] = max_mass * u.Unit(1e10*u.Msun)
    t["Max Mass Snap"] = max_mass_snapshot
    t["Current Snap Mass"] = current_mass * u.Unit(1e10*u.Msun)
    t.write(f"{groups.path_data}max_masses/{sim}_{phys}_{snap}.ecsv",
            overwrite=True)

    
