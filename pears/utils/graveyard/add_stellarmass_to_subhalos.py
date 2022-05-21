# ---------------------------------------------------------------------------- #
# the below code was used to include the median stellar mass values in the 
# abundance matching data files
# ---------------------------------------------------------------------------- #
import h5py
import numpy as np
from utils.get_groups import GetGroups
from utils.paths import SetupPaths

paths = SetupPaths()

snaps = {}

snaps["Illustris"] = np.arange(0,136,1)
snaps["TNG"] = np.arange(0,100)

for sim in ["Illustris","TNG"]:
    for snapshot in snaps[sim][::-1]:
        open_path = f"{paths.path_subhalos}{sim}_{snapshot}.hdf5"
        f = h5py.File(open_path, 'r+')

        for physics in ["dark","hydro"]:
            try:
                am_file_path = f"{sim}_{physics}_{snapshot}.hdf5"
                am = h5py.File(f"{paths.path_am_mass}{am_file_path}", "r")
                medians = np.array(am["Median Stellar Mass"])
            except OSError:
                print("AM file does not exist for", sim, snapshot, physics)
                continue

            for size in ["dwarf","massive"]:   
                try:
                    dset = f.create_dataset(f'/{physics}/{size}/Subhalo Med Stellar Mass', 
                                            shape=medians.shape,
                                            dtype=medians.dtype)
                    dset.attrs["Subhalo Med Stellar Mass"] = "Median stellar mass from abundance matching -- 1e10 Msun"
                    dset[:] = medians 

                    print("added for",sim, snapshot, physics, size)

                except RuntimeError:
                    del f[physics][size]["Subhalo Med Stellar Mass"]
                    dset = f.create_dataset(f'/{physics}/{size}/Subhalo Med Stellar Mass', 
                                            shape=medians.shape,
                                            dtype=medians.dtype)
                    dset.attrs["Subhalo Med Stellar Mass"] = "Median stellar mass from abundance matching -- 1e10 Msun"
                    dset[:] = medians 
                    print("added for",sim, snapshot, physics, size)
                
            am.close()
        f.close()
