'''
This script uses the get_groups class to save the properties of all subhalos
in the groups of interest!
'''

import h5py
import numpy as np
from utils.get_groups import GetGroups
from utils.paths import SetupPaths

paths = SetupPaths()

snaps = {}

snaps["Illustris"] = np.arange(134,136,1)
snaps["TNG"] = np.arange(99,100)

for sim in ["Illustris","TNG"]:
    for snap in snaps[sim]:
        #create the hdf5 file!
        savepath = f"{sim}_{snap}.hdf5"
        f = h5py.File(f"{paths.path_subhalos}{savepath}", 'w')

        success=False
        for phys in ["dark","hydro"]:
            am_file_path = f"{sim}_{phys}_{snap}.hdf5"
            am_masses = h5py.File(f"{paths.path_am_mass}{am_file_path}", "r")
            am_dict = {}
            for key,val in am_masses.items():
                am_dict[key] = np.array(list(val))

            for s in ["dwarf","massive"]:
                try:
                    inst = GetGroups( snapshot=snap, 
                            sim=sim,
                            physics=phys,
                            size=s)

                    maxmasses, maxmasssnaps, stellars = [], [], []

                    for i in inst.subhalo_ids:
                        ind = np.where( am_dict['Subhalo ID'] == i )[0][0]
                        maxmasses.append( am_dict['Max Mass'] )
                        maxmasssnaps.append( am_dict['Max Mass Snap'])
                        stellars.append( am_dict['Stellar Masses'] )

                    sub_dict = {"Group ID":inst.pass_numbers, 
                                "Group Mass":inst.pass_mvir, 
                                "Group Radius":inst.pass_rvir, 
                                "Nsubs":inst.pass_nsubs,
                                "Subhalo ID":inst.subhalo_ids,
                                "Subhalo Mass":inst.subhalo_masses,
                                "Subhalo Pos":inst.subhalo_pos,
                                "Subhalo Vel":inst.subhalo_vel,
                                "Subhalo Max Mass":maxmasses,
                                "Subhalo Max Mass Snap":maxmasssnaps,
                                "Subhalo Stellar Masses":stellars
                                }

                    units_dict = {
                        "Group Number":"Group Number in Subfind Catalogs", 
                        "Group Mass":"Physical mass from Group_M_TopHat200 -- 1e10 Msun", 
                        "Group Radius":"Physical radius from Group_R_TopHat200 -- kpc", 
                        "Nsubs":"Number of subhalos in group",
                        "Subhalo ID": "Subhalo ID at current snapshot",
                        "Subhalo Mass": "Subhalo mass at current snapshot -- 1e10 Msun",
                        "Subhalo Pos":"Subhalo physical position in box x,y,z -- kpc",
                        "Subhalo Vel":"Subhalo velocity in vx, vy, vz -- km/s",
                        "Subhalo Max Mass": "Maximum mass ever achieved by the subhalo -- 1e10 Msun",
                        "Subhalo Max Mass Snap": "Snapshot at which maximum halo mass is achieved",
                        "Subhalo Stellar Masses": "Stellar masses from abundance matching -- 1e10 Msun"
                        }

                    for key, val in sub_dict.items():
                        dset = f.create_dataset(f'/{phys}/{s}/{key}', 
                                                shape=val.shape,
                                                dtype=val.dtype)
                        dset.attrs[key] = units_dict[key]
                        dset[:] = val

                    success = True

                except AttributeError:
                    print(f"Cannot save {sim} {phys} {s} for snapshot {snap}")
                except OSError:
                    print(f"Cannot save {sim} {phys} {s} for snapshot {snap} - DNE")

        if success:
            #create header with simulation info
            header_dict = {"Snapshot":inst.snapshot,
                "Redshift":inst.redshift,
                "Simulation":inst.sim}

            dset = f.create_group('/Header')
            for key in header_dict.keys():
                dset.attrs[key] = header_dict[key]

            f.close()
            print(f"Saved groups at {paths.path_groups}{sim}_{snap}.hdf5")
