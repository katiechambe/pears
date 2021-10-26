import h5py
import numpy as np
from utils.get_groups import GetGroups
from utils.paths import SetupPaths

paths = SetupPaths()

snaps = {}

snaps["Illustris"] = np.arange(0,136,1)
snaps["TNG"] = np.arange(0,100)


# snaps["Illustris"] = np.arange(130,136,1)
# snaps["TNG"] = np.arange(95,100,1)

for sim in ["Illustris","TNG"]:
    for snap in snaps[sim]:

        #create the hdf5 file!
        savepath = f"{sim}_{snap}.hdf5"
        f = h5py.File(f"{paths.path_groups}{savepath}", 'w')

        success=False
        for phys in ["dark","hydro"]:
            for s in ["dwarf","massive"]:
                try:
                    inst = GetGroups( snapshot=snap, 
                            sim=sim,
                            physics=phys,
                            size=s)

                    group_dict = {"Group Number":inst.pass_numbers, 
                                "Group Mass":inst.pass_mvir, 
                                "Group Radius":inst.pass_rvir, 
                                "Nsubs":inst.pass_nsubs}

                    units_dict = {"Group Number":"Group Number in Subfind Catalogs", 
                        "Group Mass":"Physical mass from Group_M_TopHat200 -- 1e10 Msun", 
                        "Group Radius":"Physical radius from Group_R_TopHat200 -- kpc", 
                        "Nsubs":"Number of subhalos in group"}

                    for key, val in group_dict.items():
                        dset = f.create_dataset(f'/{phys}/{s}/{key}', 
                                                shape=val.shape,
                                                dtype=val.dtype)
                        dset.attrs[key] = units_dict[key]
                        dset[:] = val

                    success = True

                except AttributeError:
                    print(f"Cannot save {sim} {phys} {s} for snapshot {snap}")
                except OSError:
                    print(f"Cannot save {sim} {phys} {s} for snapshot {snap} - DNE")

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
