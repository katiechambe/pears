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
f = h5py.File(f"{paths.path_pairs}{savepath}", 'w')

redshift = subhalo_data["Header"].attrs["Redshift"]
scale = 1 / (1 + redshift)

units_dict = {
    "Group ID":"Group Number in Subfind Catalogs", 
    "Group Mass":"Physical mass from Group_M_TopHat200 -- 1e10 Msun",
    "Group Radius":"Physical radius from Group_R_TopHat200 -- kpc",
    "Group Nsubs":"Number of subhalos in group",
    "Sub1 ID":  "Subhalo ID at current snapshot",
    "Sub2 ID":  "Subhalo ID at current snapshot",
    "Sub1 Mass": "Subhalo mass at current snapshot -- 1e10 Msun",
    "Sub2 Mass": "Subhalo mass at current snapshot -- 1e10 Msun",
    "Sub1 Stellar Mass": "Stellar masses from abundance matching -- 1e10 Msun",
    "Sub2 Stellar Mass": "Stellar masses from abundance matching -- 1e10 Msun",
    "Sub1 Pos": "Subhalo physical position in box x,y,z -- kpc",
    "Sub2 Pos": "Subhalo physical position in box x,y,z -- kpc",
    "Sub1 Vel": "Subhalo velocity in vx, vy, vz -- km/s",
    "Sub2 Vel": "Subhalo velocity in vx, vy, vz -- km/s",
    "Sub1 MassType": "Mass of bound particles - gas, DM, empty, tracers, stars, BHs -- in 1e10 Msun",
    "Sub2 MassType": "Mass of bound particles - gas, DM, empty, tracers, stars, BHs -- in 1e10 Msun",
    "Separation": "Physical separation between primary and secondary in kpc",
    "RelVel": "Relative velocity between primary and secondary in km/s",
    "Stellar Mass Ratio": "Stellar mass ratio of secondary over primary", # in this case- it's by primary subhalo mass
    "Realization": "Stellar mass realization (0-1000)",
    "Sub1 BHMass": "Sum of the masses of all blackholes -- 1e10 Msun",
    "Sub2 BHMass": "Sum of the masses of all blackholes -- 1e10 Msun",
    "Sub1 BHMdot": "Instantaneous accretion rates of all blackholes -- 1e10 Msun / 0.978Gyr",
    "Sub2 BHMdot": "Instantaneous accretion rates of all blackholes -- 1e10 Msun / 0.978Gyr",
    "Sub1 SFR": "Sum of the individual SFRs of all gas cells in subhalo -- Msun / yr",
    "Sub2 SFR": "Sum of the individual SFRs of all gas cells in subhalo -- Msun / yr",
    "Sub1 SFRinRad": "Sum of SFRs of all gas cells within twice the stellar half mass radius -- Msun / yr",
    "Sub2 SFRinRad": "Sum of SFRs of all gas cells within twice the stellar half mass radius -- Msun / yr",
    "Sub1 GasMetallicity": "Mz/Mtot, where Z = any element above He within 2x stellar half mass radius -- unitless",
    "Sub2 GasMetallicity": "Mz/Mtot, where Z = any element above He within 2x stellar half mass radius -- unitless",
    "TripleFlag": "0 if no tertiary with mass ratio > 1:10 of secondary, 1 if large tertiary, 2 if other problem occurred"
    }

#create header with simulation info
header_dict = {"Snapshot":snapshot,
        "Redshift":redshift,
        "Simulation":sim}

dset = f.create_group('/Header')
for key in header_dict.keys():
    dset.attrs[key] = header_dict[key]

f.close()

class SkipFragmentedSubhalo(Exception): pass

def testthird(shortlist, stells, stel2, real):
    stells_sans_prim = np.where(stells==np.max(stells), 0, stells)
    stells_sans_two = np.where(stells_sans_prim==np.max(stells_sans_prim), 0, stells_sans_prim)
    
    if all(stells_sans_two == 0):
        print("Error in testthird")
        raise SkipFragmentedSubhalo()

    else:
        tertiary_loc = np.where(meds==np.max(stells_sans_two))
        
        if real == -1:
            stel3 = shortlist['Subhalo Med Stellar Mass'][tertiary_loc][0]
            
        else:
            stel3 = shortlist['Subhalo Stellar Masses'][tertiary_loc][0][real]

        ratio = stel3 / stel2 

        if ratio > 1/10:
            return 1
        else:
            return 2


## ensure that the dark and/or hydro file exists before operating
for phys in ["dark","hydro"]:
    try: 
        subhalo_data[phys]

    except KeyError:
        print(f"{phys} does not exist in this file: {subhalo_path}")
        break

    ## need to find pairs for dwarfs and massive

    subhalo_dict = {}
    pair_data = {}
    pair_data = {"Group ID": [],
                 "Group Mass": [],
                 "Group Radius": [],
                 "Group Nsubs": [],
                 "Sub1 ID": [],
                 "Sub2 ID": [],
                 "Sub1 Mass": [],
                 "Sub2 Mass": [],
                 "Sub1 Stellar Mass": [],
                 "Sub2 Stellar Mass": [],
                 "Sub1 Pos": [],
                 "Sub2 Pos": [],
                 "Sub1 Vel": [],
                 "Sub2 Vel": [],
                 "Separation": [],
                 "RelVel": [],
                 "Stellar Mass Ratio": [],
                 "Realization": [],
                 "TripleFlag":[],
                 "Sub1 MassType": [],
                 "Sub2 MassType": [],
                 "Sub1 BHMass": [],
                 "Sub2 BHMass": [],
                 "Sub1 BHMdot": [],
                 "Sub2 BHMdot": [],
                 "Sub1 SFR": [],
                 "Sub2 SFR": [],
                 "Sub1 SFRinRad": [],
                 "Sub2 SFRinRad": [],
                 "Sub1 GasMetallicity": [],
                 "Sub2 GasMetallicity": []
                 }

    ## create dictionary of subhalo info
    for key,val in subhalo_data[phys].items():
        subhalo_dict[key] = np.array(val)
                
    uniqueGroups = np.unique(subhalo_dict['Group ID'])

    for groupNum in uniqueGroups:
        mask = subhalo_dict['Group ID'] == groupNum
        numPassingSubs = len(subhalo_dict['Nsubs'][mask])

        try:

            if numPassingSubs >=2:
                for realization in realNumbs:
                    groupmass = subhalo_dict['Group Mass'][mask][0]
                    groupradius = subhalo_dict['Group Radius'][mask][0]
                    
                    # initializing for potential triple
                    tripleflag = 0

                    shortlist = {} # dictionary  of subhalos in group
                    for key, val in subhalo_dict.items():
                        shortlist[key] = np.array(val[mask])

                    if realization == -1:
                        meds = shortlist["Subhalo Med Stellar Mass"]
                        primary_loc = np.where(meds==np.max(meds))

                        meds_sansmax = np.where(meds==np.max(meds), 0, meds)
                        if all(meds_sansmax == 0):
                            raise SkipFragmentedSubhalo()

                        else:
                            secondary_loc = np.where(meds==np.max(meds_sansmax))
                            stel1 = shortlist['Subhalo Med Stellar Mass'][primary_loc][0]
                            stel2 =  shortlist['Subhalo Med Stellar Mass'][secondary_loc][0] 
                            
                            # test if there are 3 subhalos in the group
                            if numPassingSubs >=3:
                                meds_sans_two = np.where(meds_sansmax==np.max(meds_sansmax), 0, meds_sansmax)
                                if all(meds_sans_two == 0):
                                    raise SkipFragmentedSubhalo()

                                else:
                                    tertiary_loc = np.where(meds==np.max(meds_sans_two))
                                    stel3 = shortlist['Subhalo Med Stellar Mass'][secondary_loc][0] 

                                    ratio = stel3 / stel2 

                                    if ratio > 1/10:
                                        tripleflag = 1
                                        
                    else:
                        stells = shortlist["Subhalo Stellar Masses"][:,realization]
                        primary_loc = np.where(stells==np.max(stells))

                        stells_sansmax = np.where(stells==np.max(stells), 0, stells)
                        secondary_loc = np.where(stells==np.max(stells_sansmax))

                        stel1 = shortlist['Subhalo Stellar Masses'][primary_loc][0][realization]
                        stel2 =  shortlist['Subhalo Stellar Masses'][secondary_loc][0][realization]
                        
                        # add flag if there is a 3rd subhalo w high mass    
                        if numPassingSubs >=3:
                            stells_sans_two = np.where(stells_sansmax==np.max(stells_sansmax), 0, stells_sansmax)

                            tertiary_loc = np.where(stells==np.max(stells_sans_two))

                            stel3 = shortlist['Subhalo Stellar Masses'][tertiary_loc][0][realization]

                            ratio = stel3 / stel2 

                            if ratio > 1/10:
                                tripleflag = 1

                    id1 = shortlist['Subhalo ID'][primary_loc][0]
                    id2 = shortlist['Subhalo ID'][secondary_loc][0]
                    mass1 = shortlist['Subhalo Mass'][primary_loc][0]
                    mass2 = shortlist['Subhalo Mass'][secondary_loc][0]
                    
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
                    stlrt = stel2/stel1
                    
                    sub1masstype = shortlist['SubhaloMassType'][primary_loc][0]
                    sub2masstype = shortlist['SubhaloMassType'][secondary_loc][0]
                    
                    
                    if phys == "hydro":     
                        sub1bhm = shortlist['SubhaloBHMass'][primary_loc][0]
                        sub2bhm = shortlist['SubhaloBHMass'][secondary_loc][0]
                        sub1bhmdot = shortlist['SubhaloBHMdot'][primary_loc][0]
                        sub2bhmdot = shortlist['SubhaloBHMdot'][secondary_loc][0]
                        sub1sfr = shortlist['SubhaloSFR'][primary_loc][0]
                        sub2sfr = shortlist['SubhaloSFR'][secondary_loc][0]
                        sub1sfrrad = shortlist['SubhaloSFRinRad'][primary_loc][0]
                        sub2sfrrad = shortlist['SubhaloSFRinRad'][secondary_loc][0]
                        sub1gasmet = shortlist['SubhaloGasMetallicity'][primary_loc][0]
                        sub2gasmet = shortlist['SubhaloGasMetallicity'][secondary_loc][0]
                            
                    else:
                        sub1bhm, sub2bhm, sub1bhmdot, sub2bhmdot = 0, 0, 0, 0
                        sub1sfr, sub2sfr, sub1sfrrad, sub2sfrrad, sub1gasmet, sub2gasmet = 0, 0, 0, 0, 0, 0
                    
                    pairlist = np.array([groupNum, groupmass, groupradius, numPassingSubs, 
                                id1, id2, mass1, mass2, stel1, stel2, 
                                pos1, pos2, vel1, vel2, sep, relvel, stlrt,
                                realization, tripleflag, sub1masstype, sub2masstype,
                                sub1bhm, sub2bhm, sub1bhmdot, sub2bhmdot,
                                sub1sfr, sub2sfr, sub1sfrrad, sub2sfrrad, sub1gasmet, sub2gasmet])

                    for ind, key in enumerate(pair_data.keys()):
                        pair_data[key].append(pairlist[ind])

        except SkipFragmentedSubhalo:
            print(f"Skipping group {groupNum}")
            pass

    print(f"going to try and save {phys}")


    f = h5py.File(f"{paths.path_pairs}{savepath}", 'r+')

    for key, val in pair_data.items():
        val = np.array(val)
        dset = f.create_dataset(f'/{phys}/{key}', 
                                shape=val.shape,
                                dtype=val.dtype)
        dset.attrs[key] = units_dict[key]
        dset[:] = val
    f.close()
    print(f"successfully wrote {sim} {phys} to {savepath}")
    print(f"saved data for", pair_data.keys())

subhalo_data.close()

