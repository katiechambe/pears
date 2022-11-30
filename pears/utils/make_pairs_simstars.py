'''
This script uses sorts through the catalog of subhalos to find pairs for the 
stellar mass from the hydro sims
'''
import sys
import h5py
import numpy as np
from utils.get_groups import GetGroups
from utils.paths import SetupPaths
from vectorCorrection import vectorCorrection as vector

snapshot = int(sys.argv[1])
sim = str(sys.argv[2])

paths = SetupPaths()

## find subhalos for this simulation and snapshot #
subhalo_path = f"{sim}_{snapshot}.hdf5"
subhalo_data = h5py.File(f"{paths.path_subhalos}{subhalo_path}", "r")

savepath = f"{sim}_{snapshot}_simstars.hdf5"
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

phys = 'hydro'

pair_data = {}
pair_data = {"Group ID": [], "Group Mass": [], "Group Radius": [], "Group Nsubs": [],
             "Sub1 ID": [], "Sub2 ID": [], "Sub1 Mass": [], "Sub2 Mass": [],
             "Sub1 Stellar Mass": [], "Sub2 Stellar Mass": [], "Stellar Mass Ratio": [],
             "Sub1 Pos": [], "Sub2 Pos": [], "Separation": [],
             "Sub1 Vel": [], "Sub2 Vel": [], "RelVel": [],
             "Sub1 MassType": [], "Sub2 MassType": [],
             "Sub1 BHMass": [], "Sub2 BHMass": [],
             "Sub1 BHMdot": [], "Sub2 BHMdot": [],
             "Sub1 SFR": [], "Sub2 SFR": [],
             "Sub1 SFRinRad": [], "Sub2 SFRinRad": [],
             "Sub1 GasMetallicity": [], "Sub2 GasMetallicity": [],
             "Realization": [], "TripleFlag":[]}

unpaired_data = {}
unpaired_data = {"Group ID": [], "Group Mass": [], "Group Radius": [], "Group Nsubs": [], 
                 "Sub1 ID": [], "Sub1 Mass": [], "Sub1 Stellar Mass": [], "Sub1 MassType": [],
                 "Sub1 Pos": [], "Sub1 Vel": [],
                 "Sub1 BHMass": [], "Sub1 BHMdot": [], "Sub1 SFR": [], "Sub1 SFRinRad": [],
                 "Sub1 GasMetallicity": [], "Realization":[]}
    
subhalo_dictionary = {}
## create dictionary of subhalo info
for key,val in subhalo_data[phys].items():
    subhalo_dictionary[key] = np.array(val)

uniqueGroups = np.unique(subhalo_dictionary['Group ID'])

for groupNum in uniqueGroups:
    group_mask = subhalo_dictionary['Group ID'] == groupNum
    numPassingSubs = len(subhalo_dictionary['Nsubs'][group_mask])

    try:
        # single subhalo groups first!
        if numPassingSubs ==1 :
            groupmass = subhalo_dictionary['Group Mass'][group_mask][0]
            groupradius = subhalo_dictionary['Group Radius'][group_mask][0]
            groupNsubs = subhalo_dictionary['Nsubs'][group_mask][0]
            subid = subhalo_dictionary['Subhalo ID'][group_mask][0]
            submass = subhalo_dictionary['Subhalo Mass'][group_mask][0]

            substell = subhalo_dictionary['SubhaloMassType'][group_mask][0][4]
            
            subpos = subhalo_dictionary['Subhalo Pos'][group_mask][0]
            subvel = subhalo_dictionary['Subhalo Vel'][group_mask][0]
            subMassType = subhalo_dictionary['SubhaloMassType'][group_mask][0]
            subBHmass = subhalo_dictionary['SubhaloBHMass'][group_mask][0]
            subBHMdot = subhalo_dictionary['SubhaloBHMdot'][group_mask][0]
            subSFR = subhalo_dictionary['SubhaloSFR'][group_mask][0]
            subSFRinRad = subhalo_dictionary['SubhaloSFRinRad'][group_mask][0]
            subGasMet = subhalo_dictionary['SubhaloGasMetallicity'][group_mask][0]

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
                        "Realization": -1}

            for key in unpaired_data.keys():
                unpaired_data[key].append(single[key])
        
        # then try two+ subhalo groups!
        if numPassingSubs >=2:
        
            groupmass = subhalo_dictionary['Group Mass'][group_mask][0]
            groupradius = subhalo_dictionary['Group Radius'][group_mask][0]
            groupNsubs = subhalo_dictionary['Nsubs'][group_mask][0]


            # initializing for potential triple
            tripleflag = 0

            shortlist = {} # dictionary  of subhalos in group
            for key, val in subhalo_dictionary.items():
                shortlist[key] = np.array(val[group_mask])

            stells = shortlist["SubhaloMassType"][:,4]
            primary_loc = np.where(stells==np.max(stells))

            stells_sansmax = np.where(stells==np.max(stells), 0, stells)
            secondary_loc = np.where(stells==np.max(stells_sansmax))

            stel1 = shortlist["SubhaloMassType"][primary_loc][0][4]
            stel2 =  shortlist["SubhaloMassType"][secondary_loc][0][4]

            # add flag if there is a 3rd subhalo w high mass    
            if (numPassingSubs >=3) & (stel2 != 0) :
                stells_sans_two = np.where(stells_sansmax==np.max(stells_sansmax), 0, stells_sansmax)

                tertiary_loc = np.where(stells==np.max(stells_sans_two))

                stel3 = shortlist['SubhaloMassType'][tertiary_loc][0][4]

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

            pairs = {"Group ID": groupNum, "Group Mass": groupmass, "Group Radius": groupradius, "Group Nsubs": groupNsubs,
             "Sub1 ID": id1, "Sub2 ID": id2, "Sub1 Mass": mass1, "Sub2 Mass": mass2,
             "Sub1 Stellar Mass": stel1, "Sub2 Stellar Mass": stel2, "Stellar Mass Ratio": stlrt,
             "Sub1 Pos": pos1, "Sub2 Pos": pos2, "Separation": sep,
             "Sub1 Vel": vel1, "Sub2 Vel": vel2, "RelVel": relvel,
             "Sub1 MassType": sub1masstype, "Sub2 MassType": sub2masstype,
             "Sub1 BHMass": sub1bhm, "Sub2 BHMass": sub2bhm,
             "Sub1 BHMdot": sub1bhmdot, "Sub2 BHMdot": sub2bhmdot,
             "Sub1 SFR": sub1sfr, "Sub2 SFR": sub2sfr,
             "Sub1 SFRinRad": sub1sfrrad, "Sub2 SFRinRad": sub2sfrrad,
             "Sub1 GasMetallicity": sub1gasmet, "Sub2 GasMetallicity": sub2gasmet,
             "Realization": -1, "TripleFlag":tripleflag}


            for key in pair_data.keys():
                pair_data[key].append(pairs[key])

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
    dset.attrs[key] = units_dict[key]
    dset[:] = val 
print(f"successfully wrote unpaired halos for {sim} {phys} to {savepath}")
print(f"saved data for", unpaired_data.keys())

for key, val in pair_data.items():
    val = np.array(val)
    dset = f.create_dataset(f'/paired/{phys}/{key}', 
                            shape=val.shape,
                            dtype=val.dtype)
    dset.attrs[key] = units_dict[key]
    dset[:] = val 
print(f"successfully wrote paired halos for {sim} {phys} to {savepath}")
print(f"saved data for", pair_data.keys())

f.close()
subhalo_data.close()





