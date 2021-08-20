import numpy as np
import pandas as pd
import h5py

import hdf5libPy3 as hdf5lib
import readsubfHDF5Py3 as readsubfHDF5
import readtreeHDF5Py3 as readtreeHDF5

from .paths import SetupPaths

class TraceMergerTree(SetupPaths):

    def __init__(
        self, 
        snapshot,
        subfindID,
        sim = "Illustris",
        physics ="dark"
        ):
        """
        Identifies and pulls merger tree for a single subhalo

        Parameters
        ----------
        snapshot: int
            the number of the snapshot with the corresponding subhalo ID
        subfindID: int
            the ID number of the subhalo at the corresponding snapshot
        sim: str
            "Illustris" or "IllustrisTNG"
            to specify which simulation
        physics: str
            "dark" or "hydro"
            to specify which simulation
        """

        self.snapshot = snapshot
        self.subfindID = subfindID
        self.sim = sim
        self.physics = physics

        # TODO: fix this directory! 
        # should be a function of the simulation and the physics
        treeDirectory = 'Illustris-1-Dark-MergerTree'

        tree = readtreeHDF5.TreeDB(treeDirectory)
        branch = tree.get_main_branch( 
            snapnum, 
            subfind_id, 
            keysel=['SnapNum', 'SubhaloMass', 'SubhaloPos', 'SubhaloVel', 'SubhaloID', 'SubfindID']
            )

        self.branch = branch
        self.snaps = branch.SnapNum
        self.masses = branch.SubhaloMass
        self.positions = branch.SubhaloPos
        self.velocities = branch.SubhaloVel
        self.id = branch.SubhaloID
        self.subfindIDTree = branch.SubfindID
        
    @property
    def maxmass(self):
        """
        Max mass of the subhalo
        -- note: this only considers current and previous snapshots --

        Parameters:
        -----------
        None

        Outputs:
        --------
        maxmass: float
            the maximum mass previously achieved by a subhalo
        maxsnap: int
            the snapshot at which max mass occurs 
        maxredshift: float
            the maximum redshift at which max mass occurs
        """
        maxmass = np.max(self.masses)
        maxmass_mask =  max(self.masses)==self.masses
        maxsnap = self.snaps[maxmass_mask][0]
        return maxmass, maxsnap


