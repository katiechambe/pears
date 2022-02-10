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
                am_dict[key] = np.array(val)

            for s in ["dwarf","massive"]:
                try:
                    inst = GetGroups( snapshot=snap, 
                            sim=sim,
                            physics=phys,
                            size=s)

                    maxmasses, maxmasssnaps, stellars = ([] for i in range(3))
                    gid, gmvir, grvir, gnsubs  = ([] for i in range(4))

                    for i in inst.subhalo_ids:
                        # get the max mass and stellar mass data
                        ind = np.where( am_dict['Subhalo ID'] == i )[0][0]
                        maxmasses.append( am_dict['Max Mass'][ind] )
                        maxmasssnaps.append( am_dict['Max Mass Snap'][ind])
                        stellars.append( am_dict['Stellar Masses'][ind] )

                        # get the data about the group each subhalo belongs to
                        groupnum = am_dict['Group ID'][ind]
                        gid.append(groupnum)
                        gmvir.append(inst.mvirs_phys[groupnum])
                        grvir.append(inst.rvirs_phys[groupnum])
                        gnsubs.append(inst.nsubs[groupnum])

                    sub_dict = {"Group ID":np.array(gid), 
                                "Group Mass":np.array(gmvir), 
                                "Group Radius":np.array(grvir), 
                                "Nsubs":np.array(gnsubs),
                                "Subhalo ID":inst.subhalo_ids,
                                "Subhalo Mass":inst.subhalo_masses,
                                "Subhalo Pos":inst.subhalo_pos,
                                "Subhalo Vel":inst.subhalo_vel,
                                "Subhalo Max Mass":np.array(maxmasses),
                                "Subhalo Max Mass Snap":np.array(maxmasssnaps),
                                "Subhalo Stellar Masses":np.array(stellars)
                                }

                    units_dict = {
                        "Group ID":"Group Number in Subfind Catalogs", 
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
