import numpy as np
from utils.get_groups import GetGroups

snaps = {}

#snaps["Illustris"] = np.arange(0,136,1)
#snaps["TNG"] = np.arange(0,100)

snaps["Illustris"] = np.arange(130,136,1)
snaps["TNG"] = np.arange(95,100,1)

for sim in ["Illustris","TNG"]:
    for i in snaps[sim]:
        for x in ["dark","hydro"]:
            for s in ["dwarf","massive"]:
                try:
                    GetGroups( snapshot=i, 
                            sim=sim,
                            physics=x,
                            size=s).save_groups()
                except AttributeError:
                    print(f"Cannot save {sim} {x} {s} for snapshot {i}")
