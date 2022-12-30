'''
The median AM realization
'''

import sys
import h5py
import numpy as np
from utils.get_groups import GetGroups
from utils.paths import SetupPaths
from utils.vectorCorrection import vectorCorrection as vector

paths = SetupPaths()

snapshot = int(sys.argv[1])
sim = str(sys.argv[2])

## find subhalos for this simulation and snapshot ##
subhalo_path = f"{sim}_{snapshot}.hdf5"
subhalo_data = h5py.File(f"{paths.path_subhalos}{subhalo_path}", "r")
redshift = subhalo_data["Header"].attrs["Redshift"]
scale = 1 / (1 + redshift)

## create file for data collection ##
savepath = f"{paths.path_median}{sim}_{snapshot}.hdf5"
f = h5py.File(savepath, 'w')

header_dict = {"Snapshot": snapshot,
               "Redshift": redshift,
               "Scale": scale,
               "Simulation": sim}

dset = f.create_group('/Header')
for key in header_dict.keys():
    dset.attrs[key] = header_dict[key]
f.close()

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
          "Mass ratio": "Halo mass ratio of secondary over primary",
          "Stellar mass ratio": "Stellar mass ratio of secondary over primary",
          "Separation": "Physical separation between primary and secondary in kpc",
          "RelVel": "Relative velocity between primary and secondary in km/s"}

##################################
## start with unpaired subhalos ##
##################################
for phys in ["dark","hydro"]:
    try: 
        subhalo_data[phys] # ensures that the dark and/or hydro file exists before operating

    except KeyError:
        print(f"{phys} does not exist in this file: {subhalo_path}")
        break
    
    unpaired = {"Group ID": [],
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
             "Mass ratio":[],
             "Stellar mass ratio":[],
             "Separation":[],
             "RelVel":[] }

    ## create a dictionary of arrays from subhalo data
    subhalo_dictionary = {}
    for key,val in subhalo_data[phys].items():
        subhalo_dictionary[key] = np.array(val)
                
    uniqueGroups = np.unique(subhalo_dictionary['Group ID'])

    for groupNum in uniqueGroups:
        group_mask = subhalo_dictionary['Group ID'] == groupNum
        numPassingSubs = len(subhalo_dictionary['Nsubs'][group_mask])
        
        if numPassingSubs == 1:
            ## collect the single subhalo ##
            groupmass = subhalo_dictionary['Group Mass'][group_mask][0]
            groupradius = subhalo_dictionary['Group Radius'][group_mask][0]
            groupNsubs = subhalo_dictionary['Group Nsubs'][group_mask][0]
            subid = subhalo_dictionary['Subhalo ID'][group_mask][0]
            submass = subhalo_dictionary['Subhalo Mass'][group_mask][0]
            substell = subhalo_dictionary['Subhalo Med Stellar Mass'][group_mask][0]
            subpos = subhalo_dictionary['Subhalo Pos'][group_mask][0]
            subvel = subhalo_dictionary['Subhalo Vel'][group_mask][0]
            subMassType = subhalo_dictionary['Subhalo MassType'][group_mask][0]
            
            if phys == 'hydro':
                subBHmass = subhalo_dictionary['Subhalo BHMass'][group_mask][0]
                subBHMdot = subhalo_dictionary['Subhalo BHMdot'][group_mask][0]
                subSFR = subhalo_dictionary['Subhalo SFR'][group_mask][0]
                subSFRinRad = subhalo_dictionary['Subhalo SFRinRad'][group_mask][0]
                subGasMet = subhalo_dictionary['Subhalo GasMetallicity'][group_mask][0]
            else:
                subBHmass, subBHMdot, subSFR, subSFRinRad, subGasMet == 0, 0, 0, 0, 0
            
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
                      "Sub1 GasMetallicity": subGasMas,
                      "Mass ratio": 0,
                      "Stellar mass ratio": 0,
                      "Separation": -1,
                      "RelVel": -1 }
            
            for key in unpaired.keys():
                unpaired[key].append(single[key])
            
        else:
            try:
                ## i want to only pull out groups where there is a big halo and the next biggest (by mstar) is <1:10.
                for key, val in subhalo_dictionary.items():
                    shortlist[key] = np.array(val[mask])

                medians = shortlist["Subhalo Med Stellar Mass"] # "central halo"
                primary_loc = np.where(medians==np.max(medians))

                meds_sansmax = np.where(medians==np.max(medians), 0, medians)
                if all(meds_sansmax == 0):
                    raise SkipFragmentedSubhalo()

                else:
                    secondary_loc = np.where(meds==np.max(meds_sansmax))
                    stel1 = shortlist['Subhalo Med Stellar Mass'][primary_loc][0]
                    stel2 =  shortlist['Subhalo Med Stellar Mass'][secondary_loc][0] 

                    stellratio = stel2/stel1

                    if stellratio < 1/10:
                        groupmass = shortlist['Group Mass'][primary_loc][0]
                        groupradius = shortlist['Group Radius'][primary_loc][0]
                        groupNsubs = shortlist['Group Nsubs'][primary_loc][0]
                        subid = shortlist['Subhalo ID'][primary_loc][0]
                        submass = shortlist['Subhalo Mass'][primary_loc][0]
                        substell = shortlist['Subhalo Med Stellar Mass'][primary_loc][0]
                        subMassType = shortlist['Subhalo MassType'][primary_loc][0]

                        if phys == 'hydro':
                            subBHmass = shortlist['Subhalo BHMass'][primary_loc][0]
                            subBHMdot = shortlist['Subhalo BHMdot'][primary_loc][0]
                            subSFR = shortlist['Subhalo SFR'][primary_loc][0]
                            subSFRinRad = shortlist['Subhalo SFRinRad'][primary_loc][0]
                            subGasMet = shortlist['Subhalo GasMetallicity'][primary_loc][0]
                        else:
                            subBHmass, subBHMdot, subSFR, subSFRinRad, subGasMet == 0, 0, 0, 0, 0 

                        massratio = shortlist['Subhalo Mass'][secondary_loc][0]/shortlist['Subhalo Mass'][primary_loc][0]

                        pos1 = shortlist['Subhalo Pos'][primary_loc][0]
                        pos2 = shortlist['Subhalo Pos'][secondary_loc][0]
                        vel1 = shortlist['Subhalo Vel'][primary_loc][0]
                        vel2 = shortlist['Subhalo Vel'][secondary_loc][0]

                        # note: box size is physical units!
                        if sim == "Illustris":
                            boxsize = 106.5 # in Mpc!
                        elif sim == "TNG":
                            boxsize = 110.7 # in Mpc!

                        sep = np.linalg.norm( np.array(vector(pos1,pos2,boxsize*1000) ) )
                        relvel = np.linalg.norm(vel1-vel2)

                        single = {"Group ID": groupNum,
                                  "Group Mass": groupmass,
                                  "Group Radius": groupradius,
                                  "Group Nsubs": groupNsubs,
                                  "Sub1 ID": subid,
                                  "Sub1 Mass": submass,
                                  "Sub1 Stellar Mass": substell,
                                  "Sub1 Pos": pos1, # these are dif than in solo case
                                  "Sub1 Vel": vel1, # these are dif than in solo case
                                  "Sub1 MassType": subMassType,
                                  "Sub1 BHMass": subBHmass,
                                  "Sub1 BHMdot": subBHMdot,
                                  "Sub1 SFR": subSFR,
                                  "Sub1 SFRinRad": subSFRinRad,
                                  "Sub1 GasMetallicity": subGasMas,
                                  "Mass ratio": massratio,
                                  "Stellar mass ratio": stellratio,
                                  "Separation": sep,
                                  "RelVel": relvel}

                    else: 
                        continue

                for key in unpaired.keys():
                    unpaired[key].append(single[key])

            except SkipFragmentedSubhalo:
                print(f"skipping group {groupNum}")
                pass

    print(f"saving unpaired subhalos in {phys}. I found {len(unpaired[unpaired.keys()[0]])} unpaired subhalos!")

    f = h5py.File(savepath, 'r+')
    for key, val in unpaired.items():
        val = np.array(val)
        dset = f.create_dataset(f'/{phys}/unpaired/{key}', 
                                shape=val.shape,
                                dtype=val.dtype)
        dset.attrs[key] = units_dictionary[key]
        dset[:] = val
    f.close()
    
    print(f"successfully wrote {sim} {phys} to {savepath}")
    print(f"saved data for", unpaired.keys())

subhalo_data.close()

