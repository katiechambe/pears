""" 
Creates data file of maximum mass of all subhalos with mass >1e9
for both hydro and dark runs of either Illustris-1 or TNG

--------
How to run:
-----------
python make_max_mass_file.py snapshot# simulation

where snapshot can range between 0 and 136 for "Illustris"
and 0 to 99 for "TNG"
"""
__author__ = "Katie Chamberlain"
__status__ = "Beta - forever~"
__date__   = "October 2021"

import h5py
import sys
from utils.merger_trees import TraceMergerTree
from utils.get_groups import GetGroups

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

    save_path = f"{groups.path_data}max_masses/{sim}_{phys}_{snap}.hdf5"
    
    data_dict = {"Subhalo ID":subids, 
             "Max Mass":max_mass, 
             "Max Mass Snap":max_mass_snapshot, 
             "Current Snap Mass":current_mass}

    f = h5py.File(save_path, 'w')
    for key, val in data_dict.items():
        f[key] = val
    f.close()

    print(f"Saved at {groups.path_data}max_masses/{sim}_{phys}_{snap}.hdf5")
