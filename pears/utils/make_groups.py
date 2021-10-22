import numpy as np
from utils.get_groups import GetGroups

snaps = {}

snaps["Illustris"] = np.arange(0,136,1)
snaps["TNG"] = np.arange(0,100)

for sim in ["Illustris","TNG"]:
    for i in snaps[sim]:
        for x in ["dark","hydro"]:
            for s in ["dwarf","massive"]:
                try:
                    GetGroups( snapshot=i, 
                            sim="Illustris",
                            physics=x,
                            size=s).save_groups()
                except AttributeError:
                    print(f"Cannot save {sim} {x} {s} for snapshot {i}")
