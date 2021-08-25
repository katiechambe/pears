from utils import readtreeHDF5Py3 as readTree
from utils.paths import SetupPaths

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

        SetupPaths.__init__(self)

        self.snapshot = snapshot
        self.subfindID = subfindID
        self.sim = sim
        self.physics = physics

        # defining the simulation path from paths.py
        if self.sim == "Illustris":
            if self.physics == "dark":
                self.treepath = self.path_illustrisdark_trees
            elif self.physics == "hydro":
                self.treepath = self.path_illustrishydro_trees
                
        elif self.sim == "IllustrisTNG":
            if self.physics == "dark":
                self.treepath = self.path_tngdark_trees
            elif self.physics == "hydro":
                self.treepath = self.path_tnghydro_trees

        treeDirectory = self.treepath

        tree = readTree.TreeDB(treeDirectory)
        branch = tree.get_main_branch( 
            self.snapshot, 
            self.subfindID
            # keysel=['SnapNum', 'SubhaloMass', 'SubhaloPos', 'SubhaloVel', 'SubhaloID', 'SubfindID']
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
        maxmass = max(self.masses)
        maxmass_mask =  max(self.masses)==self.masses
        maxsnap = self.snaps[maxmass_mask][0]
        return maxmass, maxsnap


