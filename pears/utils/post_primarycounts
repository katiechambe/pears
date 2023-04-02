import h5py
import numpy as np
import matplotlib.pyplot as plt
from utils.paths import SetupPaths
from utils.get_summary_data import compile_summary
from utils.vectorCorrection import vectorCorrection as vector

paths = SetupPaths()




class EmptyFile(Exception): pass
errorprint = False
    
Illustris = {}
for it_phys in ['dark','hydro']:
    Illustris[it_phys] = {}
    z = []
    meds = []
    quarts = []
    for snapnum in np.arange(0,136):
        try:
            pair_path = f"Illustris_{snapnum}_10.hdf5"
            pair_data = h5py.File(f"{paths.path_pairs}{pair_path}", "r")
            realarray = np.unique(np.array(pair_data['pairs']['dark']['Realization']))
            
            if len(pair_data['pairs']['dark']['Group ID']) == 0:
                raise EmptyFile
            
            z.append(pair_data['Header'].attrs['Redshift'])
            
            num_prims_per_real = []
            for realization in realarray:
                # find num of unpaired prims:
                unpaired_stells = np.array(pair_data['unpaired']['dark']['Sub1 Stellar Mass'])
                unpaired_realization = np.array(pair_data['unpaired']['dark']['Realization'])

                mask_prims = (unpaired_stells > 0.01) & (unpaired_stells < 0.5)
                mask_real = unpaired_realization == realization

                num_unpairedprims = np.count_nonzero(mask_real & mask_prims)


                # find num of paired prims:
                paired_stells = np.array(pair_data['pairs']['dark']['Sub1 Stellar Mass'])
                paired_realization = np.array(pair_data['pairs']['dark']['Realization'])

                mask_prims = (paired_stells > 0.01) & (paired_stells < 0.5)
                mask_real = paired_realization == realization

                num_pairedprims = np.count_nonzero(mask_real & mask_prims)

                tot_numprims = num_unpairedprims + num_pairedprims

#                 print(f"total number of dwarf prims in real {realization} is: {tot_numprims}")
                num_prims_per_real.append(tot_numprims)

            meds.append( np.median(num_prims_per_real) )
            quarts.append( np.percentile(num_prims_per_real,[16,84]) )
        
        except KeyError:
            if errorprint: print(f'skipping {snapnum}')
            continue
            
        except EmptyFile:
            if errorprint: print(f"skipping {snapnum}, empty file")
            continue
            
    Illustris[it_phys] = {"z":z, "medians":meds, "quartiles":quarts}
            
    print(f"{it_phys} complete")

