'''
This script uses sorts through the catalog of subhalos to find pairs for the 
this realization of AM masses
'''
import sys
import h5py
import numpy as np
from utils.get_groups import GetGroups
from utils.paths import SetupPaths
from utils.vectorCorrection import vectorCorrection as vector

snapshot = int(sys.argv[1])
sim = str(sys.argv[2])
num_reals = int(sys.argv[3])

# snapshot = int(135)
# sim = str("Illustris")
# num_reals = int(sys.argv[4])

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
    "Separation": "Physical separation between primary and secondary in kpc",
    "RelVel": "Relative velocity between primary and secondary in km/s",
    "Stellar Mass Ratio": "Stellar mass ratio of secondary over primary", # in this case- it's by primary subhalo mass
    "Realization": "Stellar mass realization (0-1000)" }

#create header with simulation info
header_dict = {"Snapshot":snapshot,
        "Redshift":redshift,
        "Simulation":sim}

dset = f.create_group('/Header')
for key in header_dict.keys():
    dset.attrs[key] = header_dict[key]

f.close()

class SkipFragmentedSubhalo(Exception): pass


## ensure that the dark and/or hydro file exists before operating
for phys in ["dark","hydro"]:
    try: 
        subhalo_data[phys]

    except KeyError:
        print(f"{phys} does not exist in this file: {subhalo_path}")
        break

    ## need to find pairs for dwarfs and massive
    for size in ["dwarf","massive"]:
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
                     "Realization": []}

        ## create dictionary of subhalo info
        for key,val in subhalo_data[phys][size].items():
            subhalo_dict[key] = np.array(val)
                
        uniqueGroups = np.unique(subhalo_dict['Group ID'])

        for groupNum in uniqueGroups:
            mask = subhalo_dict['Group ID'] == groupNum
            numPassingSubs = len(subhalo_dict['Nsubs'][mask])
            
            try:

                if numPassingSubs >=2:
                    ### likely where i have to include the other realizations
                    realNumbs = np.arange(-1,num_reals+1,1)

                    for realization in realNumbs:
                        groupmass = subhalo_dict['Group Mass'][mask][0]
                        groupradius = subhalo_dict['Group Radius'][mask][0]

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


                        else:
                            stells = shortlist["Subhalo Stellar Masses"][:,realization]
                            primary_loc = np.where(stells==np.max(stells))

                            stells_sansmax = np.where(stells==np.max(stells), 0, stells)
                            secondary_loc = np.where(stells==np.max(stells_sansmax))

                            stel1 = shortlist['Subhalo Stellar Masses'][primary_loc][0][realization]
                            stel2 =  shortlist['Subhalo Stellar Masses'][secondary_loc][0][realization]

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

                        pairlist = np.array([groupNum, groupmass, groupradius, numPassingSubs, 
                                    id1, id2, mass1, mass2, stel1, stel2, 
                                    pos1, pos2, vel1, vel2, 
                                    sep, relvel, stlrt, realization])

                        for ind, key in enumerate(pair_data.keys()):
                            pair_data[key].append(pairlist[ind])

            except SkipFragmentedSubhalo:
                print(f"Skipping group {groupNum}")
                pass
                        
        print(f"going to try and save {size} {phys}")

        
        f = h5py.File(f"{paths.path_pairs}{savepath}", 'r+')

        for key, val in pair_data.items():
            val = np.array(val)
            dset = f.create_dataset(f'/{phys}/{size}/{key}', 
                                    shape=val.shape,
                                    dtype=val.dtype)
            dset.attrs[key] = units_dict[key]
            dset[:] = val
        f.close()
        print(f"successfully wrote {size} {phys} to {savepath}")

subhalo_data.close()

