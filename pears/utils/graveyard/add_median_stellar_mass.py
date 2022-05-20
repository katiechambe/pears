# ---------------------------------------------------------------------------- #
# the below code was used to include the median stellar mass values in the 
# abundance matching data files
# ---------------------------------------------------------------------------- #
import numpy as np
import sys
import h5py
from astropy.table import QTable
from utils.abundance_matching import AbundanceMatching
from utils.paths import SetupPaths

paths = SetupPaths()

for sim in ["Illustris","TNG"]:
    if sim=="Illustris":
        snaps=np.arange(0,136)
    elif sim=="TNG":
        snaps = np.arange(0,100)

    snapdata_path = f"{paths.path_snapdata}{sim}_snapdata.csv"
    snapdata = QTable.read(snapdata_path)
    
    for snapshot in snaps[::-1]:
        print(sim, snapshot)
        for physics in ["dark","hydro"]:
            redshift = snapdata['redshift'][snapdata['snapshot'] == snapshot][0]
            try:
                open_path = f"{paths.path_am_mass}{sim}_{physics}_{snapshot}.hdf5"
                f = h5py.File(open_path, 'r+')
                maxmasses = np.array(f['Max Mass'])
                maxmasses = np.where(maxmasses==0,1e-20,maxmasses)
                medians = AbundanceMatching(maxmasses*1e10,
                                            redshift, 
                                            samples=1).stellar_mass(med=True)
                medians = np.where(medians<1,0, medians/1e10)
                try:
                    f["Median Stellar Mass"] = medians
                    f.close() 
                except RuntimeError:
                    del f["Median Stellar Mass"]
                    f["Median Stellar Mass"] = medians
                    f.close()
            except OSError:
                print("this file doesn't exist")
                continue
            
            print("added for",sim, snapshot, physics)
