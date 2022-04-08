'''
This script uses sorts through the catalog of subhalos to find pairs for the 
this realization of AM masses
'''
'''
This script uses the get_groups class to save the properties of all subhalos
in the groups of interest!
'''
import sys
import h5py
import numpy as np
from utils.get_groups import GetGroups
from utils.paths import SetupPaths
from vectorCorrection import vectorCorrection as vector


snapshot = int(sys.argv[1])
sim = str(sys.argv[2])
physics = str(sys.argv[3])
num_reals = int(sys.argv[4])

paths = SetupPaths()

## find subhalos for this simulation and snapshot #
subhalo_path = f"{sim}_{snap}.hdf5"
subhalo_data = h5py.File(f"{paths.path_subhalos}{subhalo_path}", "r")

redshift = subhalo_data["Header"].attrs["Redshift"]
scale = 1 / (1 + redshift)

## ensure that the dark and/or hydro file exists before operating
pair_data = {}
for phys in ["dark","hydro"]:
    try: 
        subhalo_data[phys]

    except KeyError:
        print(f"{phys} does not exist in this file: {subhalo_path}"))
        break

    pair_data[phys] = {}
    ## need to find pairs for dwarfs and massive
    for size in ["dwarf","massive"]:
        subhalo_dict = {}
        pair_data[phys][size] = {}

        ## create dictionary of subhalo info
        for key,val in subhalo_data[phys][size].items():
            subhalo_dict[key] = np.array(val)
                
        uniqueGroups = np.unique(subhalo_dict['Group ID']))

        for groupNum in uniqueGroups:
            mask = subhalo_dict['Group ID'] == groupNum
            numSubs = subhalo_dict['Nsubs'][mask][0]

            if numSubs >=2:
                shortlist = {} # contains info about all subhalos in the group
                for key, val in subhalo_dict:
                    shortlist[key] = np.array(val[mask])
                
                # can change these to be by stellar mass instead! 
                #currently halo mass! 
                primary_loc = 0
                secondary_loc = 1

                GroupID, Group Mass, Group Radius,
                primary_ID, primary_mass, primary_stellar,
                secondary_ID, secondary_mass, secondary_stellar 
                stellarmassratio, separation, relative velocity

                id1 = shortlist['Subhalo ID'][primary_loc]
                id2 = shortlist['Subhalo ID'][secodary_loc]
                mass1 = shortlist['Subhalo Mass'][primary_loc]
                mass2 = shortlist['Subhalo Mass'][secondary_loc]
                stel1 = shortlist['Subhalo Stellar Masses'][primary_loc][0]
                stel2 =  shortlist['Subhalo Stellar Masses'][secondary_loc][0]
                pos1 = shortlist['Subhalo Pos'][primary_loc]
                pos2 = shortlist['Subhalo Pos'][secondary_loc]
                vel1 = shortlist['Subhalo Vel'][primary_loc]
                vel2 = shortlist['Subhalo Vel'][secondary_loc]

                sep = np.linalg.norm( np.array(vector(pos1,pos2,scale) )
                relvel = np.linalg.norm(vel1-vel2)


                shortlist['Group ID'][0], shortlist['Group Mass'][0], shortlist['Group Radius'], shortlist["Nsubs"], , , stel1, , 
                mass2,stel2, 
                
                stel2/stel1, sep, relvel
                
                
                
                GroupID, Group Mass, Group Radius,
                primary_ID, primary_mass, primary_stellar,
                secondary_ID, secondary_mass, secondary_stellar,
                stellarmassratio, separation, relative velocity
                

            if numSubs == 1:



            



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
        
        else:
            print("Closing file")
            f.close()
            print(f"Something went wrong at {paths.path_groups}{sim}_{snap}.hdf5")



for x in ['test1','test2']:
    for phys in ['hydro','dark']:
        l=0
        try:
            subhalo_data[phys]
            print('here')
            l+=1
        except KeyError:
            print('DNE')
            continue
        l += 1
        print(l)
    print(x)




subhalo_data.close()
savepath = f"{sim}_{snap}.hdf5"
f = h5py.File(f"{paths.path_pairs}{savepath}", 'a')
f.close()
