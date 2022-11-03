__author__ = "Katie Chamberlain"
__status__ = "Beta - forever~"
__date__   = "September 2021"

import numpy as np
from utils.read_group_cats import ReadCats
from utils.merger_trees import TraceMergerTree
from utils.get_groups import GetGroups
from utils.paths import SetupPaths

import sys
# from astropy.table import QTable
# import astropy.units as u

snap = int(sys.argv[1])
simulation = str(sys.argv[2])

print(snap)
print(simulation)
print(type(snap))
print(type(simulation))

print(simulation=="Illustris")
