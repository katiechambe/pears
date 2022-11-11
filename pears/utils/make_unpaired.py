'''
This script uses sorts through the catalog of subhalos to find pairs for the 
this realization of AM masses
'''
import sys
import h5py
import numpy as np
from utils.get_groups import GetGroups
from utils.paths import SetupPaths
from vectorCorrection import vectorCorrection as vector

snapshot = int(sys.argv[1])
sim = str(sys.argv[2])
num_reals = int(sys.argv[3])

realNumbs = np.arange(-1,num_reals+1,1)
paths = SetupPaths()

## find subhalos for this simulation and snapshot #
subhalo_path = f"{sim}_{snapshot}.hdf5"
subhalo_data = h5py.File(f"{paths.path_subhalos}{subhalo_path}", "r")

savepath = f"{sim}_{snapshot}_10.hdf5"
# f = h5py.File(f"{paths.path_pairs}{savepath}", 'w')

redshift = subhalo_data["Header"].attrs["Redshift"]
scale = 1 / (1 + redshift)

units_dictionary = {"Group ID": "Group Number in Subfind Catalogs", 
          "Group Mass": "Physical mass from Group_M_TopHat200 -- 1e10 Msun",
          "Group Radius": "Physical radius from Group_R_TopHat200 -- kpc",
          "Group Nsubs": "Number of subhalos in group",
          "Sub1 ID": "Subhalo ID at current snapshot",
          "Sub1 Mass": "Subhalo mass at current snapshot -- 1e10 Msun",
          "Sub1 Stellar Mass": "Stellar masses from abundance matching -- 1e10 Msun",
          "Sub1 Pos": "Subhalo physical position in box x,y,z -- kpc",
          "Sub1 Vel": "Subhalo velocity in vx, vy, vz -- km/s",
          "Sub1 MassType": "Mass of bound particles - gas, DM, empty, tracers, stars, BHs -- in 1e10 Msun",
          "Sub1 BHMass": "Sum of the masses of all blackholes -- 1e10 Msun",
          "Sub1 BHMdot": "Instantaneous accretion rates of all blackholes -- 1e10 Msun / 0.978Gyr",
          "Sub1 SFR": "Sum of the individual SFRs of all gas cells in subhalo -- Msun / yr",
          "Sub1 SFRinRad": "Sum of SFRs of all gas cells within twice the stellar half mass radius -- Msun / yr",
          "Sub1 GasMetallicity": "Mz/Mtot, where Z = any element above He within 2x stellar half mass radius -- unitless",
          "Realization": "Stellar mass realization (0-1000)"
}

# header already exists since running this file after the pair creation file

class SkipFragmentedSubhalo(Exception): pass

## ensure that the dark and/or hydro file exists before operating
for phys in ["dark","hydro"]:
    try: 
        subhalo_data[phys]

    except KeyError:
        print(f"{phys} does not exist in this file: {subhalo_path}")
        break

    ## need to find pairs for dwarfs and massive
    unpaired_data = {}
    unpaired_data = {"Group ID": [],
                 "Group Mass": [],
                 "Group Radius": [],
                 "Group Nsubs": [],
                 "Sub1 ID": [],
                 "Sub1 Mass": [],
                 "Sub1 Stellar Mass": [],
                 "Sub1 Pos": [],
                 "Sub1 Vel": [],
                 "Sub1 MassType": [],
                 "Sub1 BHMass": [],
                 "Sub1 BHMdot": [],
                 "Sub1 SFR": [],
                 "Sub1 SFRinRad": [],
                 "Sub1 GasMetallicity": [],
                 "Realization":[]}

    ## create a dictionary of arrays from subhalo data
    subhalo_dictionary = {}
    for key,val in subhalo_data[phys].items():
        subhalo_dictionary[key] = np.array(val)
                
    uniqueGroups = np.unique(subhalo_dictionary['Group ID'])
    
    for groupNum in uniqueGroups:
        group_mask = subhalo_dictionary['Group ID'] == groupNum
        numPassingSubs = len(subhalo_dictionary['Nsubs'][group_mask])

        try:

            if numPassingSubs ==1 :
                for realization in realNumbs:
                    groupmass = subhalo_dictionary['Group Mass'][group_mask][0]
                    groupradius = subhalo_dictionary['Group Radius'][group_mask][0]
                    groupNsubs = subhalo_dictionary['Nsubs'][group_mask][0]
                    subid = subhalo_dictionary['Subhalo ID'][group_mask][0]
                    submass = subhalo_dictionary['Subhalo Mass'][group_mask][0]
                    
                    if realization == -1: 
                        substell = subhalo_dictionary['Subhalo Med Stellar Mass'][group_mask][0]
                    else:
                        substell = subhalo_dictionary['Subhalo Stellar Masses'][group_mask][0][realization]
                    
                    subpos = subhalo_dictionary['Subhalo Pos'][group_mask][0]
                    subvel = subhalo_dictionary['Subhalo Vel'][group_mask][0]
                    subMassType = subhalo_dictionary['SubhaloMassType'][group_mask][0]

                    if phys == 'hydro':
                        subBHmass = subhalo_dictionary['SubhaloBHMass'][group_mask][0]
                        subBHMdot = subhalo_dictionary['SubhaloBHMdot'][group_mask][0]
                        subSFR = subhalo_dictionary['SubhaloSFR'][group_mask][0]
                        subSFRinRad = subhalo_dictionary['SubhaloSFRinRad'][group_mask][0]
                        subGasMet = subhalo_dictionary['SubhaloGasMetallicity'][group_mask][0]
                    else:
                        subBHmass, subBHMdot, subSFR, subSFRinRad, subGasMet = 0, 0, 0, 0, 0

                    single = {"Group ID": groupNum,
                              "Group Mass": groupmass,
                              "Group Radius": groupradius,
                              "Group Nsubs": groupNsubs,
                              "Sub1 ID": subid,
                              "Sub1 Mass": submass,
                              "Sub1 Stellar Mass": substell,
                              "Sub1 Pos": subpos,
                              "Sub1 Vel": subvel,
                              "Sub1 MassType": subMassType,
                              "Sub1 BHMass": subBHmass,
                              "Sub1 BHMdot": subBHMdot,
                              "Sub1 SFR": subSFR,
                              "Sub1 SFRinRad": subSFRinRad,
                              "Sub1 GasMetallicity": subGasMet,
                              "Realization": realization}
                                     
                    for key in unpaired_data.keys():
                        unpaired_data[key].append(single[key])

        except SkipFragmentedSubhalo:
            print(f"Skipping group {groupNum}")
            pass

    print(f"going to try and save {phys}")


    f = h5py.File(f"{paths.path_pairs}{savepath}", 'r+')
    
    for key, val in unpaired_data.items():
        val = np.array(val)
        dset = f.create_dataset(f'/unpaired/{phys}/{key}', 
                                shape=val.shape,
                                dtype=val.dtype)
        dset.attrs[key] = units_dictionary[key]
        dset[:] = val 
    
    f.close()
    print(f"successfully wrote unpaired halos for {sim} {phys} to {savepath}")
    print(f"saved data for", unpaired_data.keys())

subhalo_data.close()

