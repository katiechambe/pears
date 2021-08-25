""" 
This script reads the subhalo from catalogs

Dependencies:
-------------
ReadCats class from read_group_cats.py

Note:
-----
outputs from this script are in SIMULATION units!

"""

__author__ = "Katie Chamberlain"
__status__ = "Beta -- forever."

from utils.paths import SetupPaths
import utils.readsubfHDF5Py3 as readSub

class ReadCats:

    def __init__(
        self,
        snapshot,
        sim = "Illustris",
        physics ="dark"
        ):
        """
        Reads group catalog for correct simulation and physics

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
        kwargs: dict
            little_h: float
                definition of little h, default h=0.702

        Outputs:
        --------
        outputs from this script are in SIMULATION units!
        """
        # TODO: get this to actually run on IllustrisTNG lmao kill me

        self.snapshot = snapshot 
        self.sim = sim
        self.physics = physics

        SetupPaths.__init__(self)
        
        # defining the simulation path from paths.py
        if self.sim == "Illustris":
            if self.physics == "dark":
                self.catpath = self.path_illustrisdark
            elif self.physics == "hydro":
                self.catpath = self.path_illustrishydro
                
        elif self.sim == "TNG":
            if self.physics == "dark":
                self.catpath = self.path_tngdark
            elif self.physics == "hydro":
                self.catpath = self.path_tnghydro

        keysel = [
            'GroupPos', 'Group_M_TopHat200', 
            'Group_R_TopHat200', 'GroupNsubs'
            ]

        self.catalog = readSub.subfind_catalog(
            basedir=self.catpath,
            snapnum=self.snapshot,
            keysel=keysel
            )

        self.groupPos = self.catalog.GroupPos
        self.mvirs = self.catalog.Group_M_TopHat200
        self.rvirs = self.catalog.Group_R_TopHat200
        self.nsubs = self.catalog.GroupNsubs
