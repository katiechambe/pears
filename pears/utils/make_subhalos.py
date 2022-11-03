'''
This script uses the get_groups class to save the properties of all subhalos
in the groups of interest!
'''

import h5py
import numpy as np
import sys
from utils.get_groups import GetGroups
from utils.paths import SetupPaths

paths = SetupPaths()

# snaps = {}

snap = int(sys.argv[1])
sim = str(sys.argv[2])

#create the hdf5 file!
savepath = f"{sim}_{snap}.hdf5"
f = h5py.File(f"{paths.path_subhalos}{savepath}", 'w')

# list of all kwargs, note that dark sims will not have last 5 fields
keysel_list = ['GroupPos','Group_M_TopHat200', 
            'Group_R_TopHat200','GroupNsubs',
            'GroupFirstSub',
            'SubhaloGrNr','SubhaloMass',
            'SubhaloPos','SubhaloVel' ,
            'SubhaloBHMass','SubhaloBHMdot',
            'SubhaloSFR','SubhaloSFRinRad',
            'SubhaloGasMetallicity','SubhaloMassType']

kwargs = {"group_mass_min": 8,
          "group_mass_max": 650,
          "keysel":keysel_list} 

base_units = {
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
            "Subhalo Stellar Masses": "Stellar masses from abundance matching -- 1e10 Msun",
            "Subhalo Med Stellar Mass":"Median stellar mass from abundance matching -- 1e10 Msun",
            "SubhaloMassType": "Mass of bound particles - gas, DM, empty, tracers, stars, BHs -- in 1e10 Msun"
            }

hydro_units = {
            "SubhaloBHMass": "Sum of the masses of all blackholes in this subhalo -- 1e10 Msun",
            "SubhaloBHMdot": "Instantaneous accretion rates of all blackholes in subhalo -- 1e10 Msun / 0.978Gyr",
            "SubhaloSFR": "Sum of the individual SFRs of all gas cells in subhalo -- Msun / yr",
            "SubhaloSFRinRad": "Sum of SFRs of all gas cells within twice the stellar half mass radius -- Msun / yr",
            "SubhaloGasMetallicity": "Mz/Mtot, where Z = any element above He within 2x stellar half mass radius -- unitless"
            }

success=False
for phys in ["dark","hydro"]:
    try:
        am_file_path = f"{sim}_{phys}_{snap}.hdf5"
        am_masses = h5py.File(f"{paths.path_am_mass}{am_file_path}", "r")
        am_dict = {}
        for key,val in am_masses.items():
            am_dict[key] = np.array(val)

        try:
            inst = GetGroups( snapshot=snap, 
                    sim=sim,
                    physics=phys,
                    **kwargs)

            maxmasses, maxmasssnaps = ([] for i in range(2))
            medstellars, stellars =[],[]
            gid, gmvir, grvir, gnsubs  = ([] for i in range(4))

            for i in inst.subhalo_ids:
                # get the max mass and stellar mass data
                ind = np.where( am_dict['Subhalo ID'] == i )[0][0]
                maxmasses.append( am_dict['Max Mass'][ind] )
                maxmasssnaps.append( am_dict['Max Mass Snap'][ind])
                stellars.append( am_dict['Stellar Masses'][ind] )
                medstellars.append( am_dict['Median Stellar Mass'][ind] )

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
                    "Subhalo Stellar Masses":np.array(stellars),
                    "Subhalo Med Stellar Mass":np.array(medstellars),
                    "SubhaloMassType":np.array(inst.catalog.SubhaloMassType[inst.subhalo_ids])/inst.little_h,
                    }
            units_dict = base_units

            if phys == "hydro":
                sub_dict = { **sub_dict,
                            "SubhaloBHMass": np.array( inst.catalog.SubhaloBHMass[inst.subhalo_ids]/inst.little_h),
                            "SubhaloBHMdot": np.array( inst.catalog.SubhaloBHMdot[inst.subhalo_ids] ),
                            "SubhaloSFR": np.array( inst.catalog.SubhaloSFR[inst.subhalo_ids] ),
                            "SubhaloSFRinRad": np.array( inst.catalog.SubhaloSFRinRad[inst.subhalo_ids] ),
                            "SubhaloGasMetallicity":np.array( inst.catalog.SubhaloGasMetallicity[inst.subhalo_ids] )
                           }

                units_dict = { **base_units, **hydro_units }

            for key, val in sub_dict.items():
                dset = f.create_dataset(f'/{phys}/{key}', 
                                        shape=val.shape,
                                        dtype=val.dtype)
                dset.attrs[key] = units_dict[key]
                dset[:] = val

            success = True


        except AttributeError:
            print(f"Cannot save {sim} {phys} for snapshot {snap}")
        except OSError:
            print(f"Cannot save {sim} {phys} for snapshot {snap} - DNE")
    except OSError:
        print(f"Cannot find AM file for {sim} {phys} for snapshot {snap}")
        print("Closed file")

if success:
    #create header with simulation info
    header_dict = {"Snapshot":inst.snapshot,
        "Redshift":inst.redshift,
        "Simulation":inst.sim}

    dset = f.create_group('/Header')
    for key in header_dict.keys():
        dset.attrs[key] = header_dict[key]

    f.close()
    print(f"Saved subhalos at {paths.path_subhalos}{sim}_{snap}.hdf5")

else:
    print("Closing file")
    f.close()
    print(f"Something went wrong at {paths.path_subhalos}{sim}_{snap}.hdf5")
