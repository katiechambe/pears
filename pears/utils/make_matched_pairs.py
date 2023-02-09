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
paths = SetupPaths()

snapshot = int(sys.argv[1])
sim = str(sys.argv[2])
num_reals = int(sys.argv[3])

## get pairs
pair_path = f"{sim}_{snapshot}_{num_reals}.hdf5"
pair_data = h5py.File(f"{paths.path_pairs}{pair_path}", "r")

#################################################
# defining save location & creating header data #
#################################################
savepath = f"{sim}_{snapshot}_{num_reals}_matched.hdf5"
f = h5py.File(f"{paths.path_data}matched/Vicente/{savepath}", 'w')
redshift = pair_data["Header"].attrs["Redshift"]
scale = 1 / (1 + redshift)
header_dict = {"Snapshot":snapshot, "Redshift":redshift, "Simulation":sim}
dset = f.create_group('/Header')
for key in header_dict.keys():
    dset.attrs[key] = header_dict[key]
f.close()

## define units of saved quantities
units_dictionary = {
          ### Group info
        "Group ID":"Group Number in Subfind Catalogs", 
        "Group Mass":"Physical mass from Group_M_TopHat200 -- 1e10 Msun",
        "Group Radius":"Physical radius from Group_R_TopHat200 -- kpc",
        "Group Nsubs":"Number of subhalos in group",
          ### Subhalo info
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
          ### Relative vals
        "Separation": "Physical separation between primary and secondary in kpc",
        "Comoving Separation":"Comoving separation between primary and secondary in ckpc",
        "RelVel": "Relative velocity between primary and secondary in km/s",
        "Stellar Mass Ratio": "Stellar mass ratio of secondary over primary", # in this case- it's by primary subhalo mass
        "Realization": "Stellar mass realization (0-1000)",
        "TripleFlag": "0 if no tertiary with mass ratio > 1:10 of secondary, 1 if large tertiary, 2 if other problem occurred"
        }



#########################
# open matched catalogs #
#########################
# reading in Nelson matched catalog
# match_Npath = h5py.File(f"{paths.path_tngmatch_N}", "r")
# matchedN = match_Npath[f'Snapshot_{snapshot}']
# darkidsN = np.array(matchedN['SubhaloIndexDark_SubLink'])

## Vicente's matched catalog
matchedV = h5py.File(f"{paths.path_tngmatch_V}dm_match_{str(snapshot).zfill(3)}.hdf5", "r")
darkidsV = np.array(matchedV["IndexTo"])
hydroidsV = np.array(matchedV["IndexFrom"])

################################################################
# open and create arrays of the subhalo data for fast indexing #
################################################################
subhalo_path = f"{sim}_{snapshot}.hdf5"
subhalo_data = h5py.File(f"{paths.path_subhalos}{subhalo_path}", "r")

dark_subhalo_dict = {}
for key,val in subhalo_data["dark"].items():
    dark_subhalo_dict[key] = np.array(val)


##################################
# start with the unpaired halos! #
##################################
unpaired = pair_data['unpaired']["hydro"]

unpair_hydro = {}
unpair_hydro = {"Group ID": [],
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

unpair_dark_V = {}
for key in unpair_hydro.keys():
    unpair_dark_V[key] = []


## create a dictionary of arrays from subhalo data
unpair_dictionary = {}
for key,val in unpaired.items():
    unpair_dictionary[key] = np.array(val)
    
missingmatch = []   
missingdark = []
for index in range(len(unpair_dictionary['Sub1 ID'])):
    if index%1000 ==0:
        print(f"unpair count: {index} of {len(unpair_dictionary['Sub1 ID'])}")
    
    id_hydro = unpair_dictionary['Sub1 ID'][index]

    hydro_single = {}
    for key in unpair_hydro.keys():
        hydro_single[key] = unpair_dictionary[key][index]

    real = hydro_single["Realization"]
    
        ## check if the hydro subhalo has a match in dark
    if id_hydro not in hydroidsV:
#         print(f"{id_hydro} not in Vicente catalog")
        missingmatch.append(id_hydro)
        continue
        
    id_dark_V = darkidsV[np.where(hydroidsV == id_hydro)[0][0]]
    
    if id_dark_V not in dark_subhalo_dict["Subhalo ID"]:
#         print(f"{id_dark_V} is not in dark catalog")
        missingdark.append(id_dark_V)
        continue
        
    dark_single_V = {}
    index_V = np.where(dark_subhalo_dict['Subhalo ID'] == id_dark_V)[0][0]
    
    if real == -1:
        dark_V_stell = dark_subhalo_dict['Subhalo Med Stellar Mass'][index_V]
    else:
        dark_V_stell = dark_subhalo_dict['Subhalo Stellar Masses'][index_V][real]
    
    dark_single_V = {"Group ID": dark_subhalo_dict['Group ID'][index_V],
                     "Group Mass": dark_subhalo_dict['Group Mass'][index_V],
                     "Group Radius": dark_subhalo_dict["Group Radius"][index_V],
                     "Group Nsubs": dark_subhalo_dict['Nsubs'][index_V],
                     "Sub1 ID": dark_subhalo_dict['Subhalo ID'][index_V],
                     "Sub1 Mass": dark_subhalo_dict['Subhalo Mass'][index_V],
                     "Sub1 Stellar Mass": dark_V_stell,
                     "Sub1 Pos": dark_subhalo_dict['Subhalo Pos'][index_V],
                     "Sub1 Vel": dark_subhalo_dict['Subhalo Vel'][index_V],
                     "Sub1 MassType": dark_subhalo_dict['SubhaloMassType'][index_V],
                     "Sub1 BHMass": 0,
                     "Sub1 BHMdot": 0,
                     "Sub1 SFR": 0,
                     "Sub1 SFRinRad": 0,
                     "Sub1 GasMetallicity": 0,
                     "Realization":real}
    
    for key in unpair_hydro.keys():
        unpair_hydro[key].append(hydro_single[key])
        unpair_dark_V[key].append(dark_single_V[key])
    
    
print(f"unpaired data structures built. saving snapshot:{snapshot}.")

f = h5py.File(f"{paths.path_data}matched/Vicente/{savepath}", 'r+')

for key, val in unpair_hydro.items():
    val = np.array(val)
    dset = f.create_dataset(f'/unpaired/hydro/{key}', 
                            shape=val.shape,
                            dtype=val.dtype)
    dset.attrs[key] = units_dictionary[key]
    dset[:] = val 
    
for key, val in unpair_dark_V.items():
    val = np.array(val)
    dset = f.create_dataset(f'/unpaired/dark match/{key}', 
                            shape=val.shape,
                            dtype=val.dtype)
    dset.attrs[key] = units_dictionary[key]
    dset[:] = val 

f.close()
print(f"successfully wrote unpaired halos for {sim} {snapshot} to {savepath}")
print(f"total number of hydro subhalos: {len(np.array(unpair_dictionary['Sub1 ID']))}")
print(f"total number of unique hydro subhalos: {len(np.unique(np.array(unpair_dictionary['Sub1 ID'])))}")
print(f"total number of halos not in matched catalog: {len(np.unique(np.array(missingmatch)))}")
print(f"total number of halos not in dark catalog: {len(np.unique(np.array(missingdark)))}")
print(f"saved data for", unpair_dark_V.keys())

    
        
######################
# paired halos next! #
######################
paired = pair_data["pairs"]["hydro"]

pair_hydro = {}
pair_hydro = {"Group ID": [],
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
             "Sub2 GasMetallicity": [],
             "Separation": [],
             "Comoving Separation": [],
             "RelVel": [],
             "Stellar Mass Ratio": [],
             "Realization": [],
             "TripleFlag":[]}

pair_dark_V = {}
for key in pair_hydro.keys():
    pair_dark_V[key] = []

## create a dictionary of arrays from subhalo data
pair_dictionary = {}
for key,val in paired.items():
    pair_dictionary[key] = np.array(val)
    
missingmatch = []   
missingdark = []
for index in range(len(pair_dictionary['Sub1 ID'])):
    if index%10000 ==0:
        print(f"pair count: {index} of {len(pair_dictionary['Sub1 ID'])}")
    
    id_hydro = pair_dictionary['Sub1 ID'][index]
    id_hydro2 = pair_dictionary['Sub2 ID'][index]

    hydro_singlepair = {}
    for key in pair_hydro.keys():
        hydro_singlepair[key] = pair_dictionary[key][index]

    real = hydro_singlepair["Realization"]
    
    ## check if the hydro subhalo has a match in dark
    if (id_hydro not in hydroidsV) or (id_hydro2 not in hydroidsV):
      #  print(f"{id_hydro} or {id_hydro2} not in Vicente catalog")
        missingmatch.append(id_hydro2)
        continue
        
    id_dark_V = darkidsV[np.where(hydroidsV == id_hydro)[0][0]]
    id_dark_V2 = darkidsV[np.where(hydroidsV == id_hydro2)[0][0]]
    
    if (id_dark_V not in dark_subhalo_dict["Subhalo ID"]) or (id_dark_V2 not in dark_subhalo_dict["Subhalo ID"]):
#         print(f"{id_dark_V} or {id_dark_V2} is not in dark catalog")
        missingdark.append(id_dark_V2)
        continue
        
    dark_singlepair_V = {}
    index_V = np.where(dark_subhalo_dict['Subhalo ID'] == id_dark_V)[0][0]
    index_V2 = np.where(dark_subhalo_dict['Subhalo ID'] == id_dark_V2)[0][0]
    
    if real == -1:
        dark_V_stell = dark_subhalo_dict['Subhalo Med Stellar Mass'][index_V]
        dark_V_stell2 = dark_subhalo_dict['Subhalo Med Stellar Mass'][index_V2]
    else:
        dark_V_stell = dark_subhalo_dict['Subhalo Stellar Masses'][index_V][real]
        dark_V_stell2 = dark_subhalo_dict['Subhalo Stellar Masses'][index_V2][real]
      
    pos1 = dark_subhalo_dict['Subhalo Pos'][index_V]
    pos2 = dark_subhalo_dict['Subhalo Pos'][index_V2]
    copos1 = pos1/scale
    copos2 = pos2/scale

    vel1 = dark_subhalo_dict['Subhalo Vel'][index_V]
    vel2 = dark_subhalo_dict['Subhalo Vel'][index_V2]

    if sim == "Illustris":
        little_h = 0.704

    elif sim == "TNG":
        little_h = 0.6774

    boxsize = 75000.0 # in cMpc/h
    boxsize_phys = boxsize * scale / little_h    
    boxsize_co = boxsize / little_h    

    sep = np.linalg.norm( np.array(vector(pos1,pos2,boxsize_phys) ) )
    cosep = np.linalg.norm( np.array(vector(copos1,copos2,boxsize_co) ) )

    relvel = np.linalg.norm(vel1-vel2)
    stlrt = dark_V_stell2/dark_V_stell
        
    dark_singlepair_V = {"Group ID": dark_subhalo_dict['Group ID'][index_V],
             "Group Mass": dark_subhalo_dict['Group Mass'][index_V],
             "Group Radius": dark_subhalo_dict["Group Radius"][index_V],
             "Group Nsubs": dark_subhalo_dict['Nsubs'][index_V],
             "Sub1 ID": dark_subhalo_dict['Subhalo ID'][index_V],
             "Sub2 ID": dark_subhalo_dict['Subhalo ID'][index_V2],
             "Sub1 Mass": dark_subhalo_dict['Subhalo Mass'][index_V],
             "Sub2 Mass": dark_subhalo_dict['Subhalo Mass'][index_V2],
             "Sub1 Stellar Mass": dark_V_stell,
             "Sub2 Stellar Mass": dark_V_stell2,
             "Sub1 Pos": pos1,
             "Sub2 Pos": pos2,
             "Sub1 Vel": vel1,
             "Sub2 Vel": vel2,
             "Sub1 MassType": dark_subhalo_dict['SubhaloMassType'][index_V],
             "Sub2 MassType": dark_subhalo_dict['SubhaloMassType'][index_V2],
             "Sub1 BHMass": 0,
             "Sub2 BHMass": 0,
             "Sub1 BHMdot": 0,
             "Sub2 BHMdot": 0,
             "Sub1 SFR": 0,
             "Sub2 SFR": 0,
             "Sub1 SFRinRad": 0,
             "Sub2 SFRinRad": 0,
             "Sub1 GasMetallicity": 0,
             "Sub2 GasMetallicity": 0,
             "Separation": sep,
             "Comoving Separation": cosep,
             "RelVel": relvel,
             "Stellar Mass Ratio": stlrt,
             "Realization": real,
             "TripleFlag":hydro_singlepair["TripleFlag"]}
    
    for key in pair_hydro.keys():
        pair_hydro[key].append(hydro_singlepair[key])
        pair_dark_V[key].append(dark_singlepair_V[key])
    
    
print(f"pair data structures built. saving snapshot:{snapshot}.")

f = h5py.File(f"{paths.path_data}matched/Vicente/{savepath}", 'r+')

for key, val in pair_hydro.items():
    val = np.array(val)
    dset = f.create_dataset(f'/paired/hydro/{key}', 
                            shape=val.shape,
                            dtype=val.dtype)
    dset.attrs[key] = units_dictionary[key]
    dset[:] = val 
    
for key, val in pair_dark_V.items():
    val = np.array(val)
    dset = f.create_dataset(f'/paired/dark match/{key}', 
                            shape=val.shape,
                            dtype=val.dtype)
    dset.attrs[key] = units_dictionary[key]
    dset[:] = val 

f.close()
print(f"successfully wrote paired halos for {sim} {snapshot} to {savepath}")
print(f"total number of hydro subhalos: {len(np.array(pair_dictionary['Sub1 ID']))}")
print(f"total number of unique hydro subhalos: {len(np.unique(np.array(pair_dictionary['Sub1 ID'])))}")
print(f"total number of halos not in matched catalog: {len(np.unique(np.array(missingmatch)))}")
print(f"total number of halos not in dark catalog: {len(np.unique(np.array(missingdark)))}")
print(f"saved data for", pair_dark_V.keys())
        
        
        
        
        
        
        








# subhalodat = subhalo_data["hydro"]

# subhalo_dict = {}
# pair_data = {}
# pair_data = {"Group ID": [],
#              "Group Mass": [],
#              "Group Radius": [],
#              "Group Nsubs": [],
#              "Sub1 ID": [],
#              "Sub2 ID": [],
#              "Sub1 Mass": [],
#              "Sub2 Mass": [],
#              "Sub1 Stellar Mass": [],
#              "Sub2 Stellar Mass": [],
#              "Sub1 Pos": [],
#              "Sub2 Pos": [],
#              "Sub1 Vel": [],
#              "Sub2 Vel": [],
#              "Sub1 MassType": [],
#              "Sub2 MassType": [],
#              "Sub1 BHMass": [],
#              "Sub2 BHMass": [],
#              "Sub1 BHMdot": [],
#              "Sub2 BHMdot": [],
#              "Sub1 SFR": [],
#              "Sub2 SFR": [],
#              "Sub1 SFRinRad": [],
#              "Sub2 SFRinRad": [],
#              "Sub1 GasMetallicity": [],
#              "Sub2 GasMetallicity": []
#              "Separation": [],
#              "Comoving Separation": [],
#              "RelVel": [],
#              "Stellar Mass Ratio": [],
#              "Realization": [],
#              "TripleFlag":[]
#              }

# ## create dictionary of subhalo info
# for key,val in ssubhalodat.items():
#     subhalo_dict[key] = np.array(val)

# uniqueGroups = np.unique(subhalo_dict['Group ID'])

# for groupNum in uniqueGroups:
#         mask = subhalo_dict['Group ID'] == groupNum
#         numPassingSubs = len(subhalo_dict['Nsubs'][mask])
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        

