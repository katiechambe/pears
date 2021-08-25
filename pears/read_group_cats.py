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
        """
        # TODO: get this to actually run on IllustrisTNG lmao kill me

        self.snapshot = snapshot 
        self.h = 0.704 # little h is not a fxn of z, by definition

        SetupPaths.__init__(self)
        
        self.sim = sim
        self.physics = physics

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

        keysel = ['GroupFirstSub','GroupPos','SubhaloMass','SubhaloPos',
                  'SubhaloGrNr', 'Group_M_TopHat200','Group_R_TopHat200',
                  'SubhaloVel','GroupNsubs']

        self.catalog = readSub.subfind_catalog(
            basedir=self.catpath,
            snapnum=self.snapshot,
            keysel=keysel
            )

        self.inds = self.catalog.GroupFirstSub
        self.groupPos = self.catalog.GroupPos
        self.mvirs = self.catalog.Group_M_TopHat200
        self.rvirs = self.catalog.Group_R_TopHat200
        self.nsubs = self.catalog.GroupNsubs
        self.submass = self.catalog.SubhaloMass
        self.subpos = self.catalog.SubhaloPos
        self.subvel = self.catalog.SubhaloVel
        self.subgr = self.catalog.SubhaloGrNr

        # TODO: must define scale! figure out how to read in the redshift!

        self.mvirs = mvirs/h # to get mass in Msun
        self.submass = submass/h # to get mass in Msun
        self.rvirs = rvirs * scale # to get rvir in kpc physical


        
##############################################
# defining each column of the data structure #
##############################################
inds = obj.GroupFirstSub
groupPos = obj.GroupPos
mvirs = obj.Group_M_TopHat200
mvirs = mvirs/h # to get mass in Msun
rvirs = obj.Group_R_TopHat200
rvirs = rvirs * scale # to get rvir in kpc physical
nsubs = obj.GroupNsubs

submass = obj.SubhaloMass
submass = submass/h # to get mass in Msun
subpos = obj.SubhaloPos
subvel = obj.SubhaloVel
subgr = obj.SubhaloGrNr


###################################
#####     Group Mass Cuts     #####
###################################
lowerGroupMass = 8 # 8e10 Msun
upperGroupMass = 50 # 5e11 Msun
groupMask = np.where((mvirs > lowerGroupMass) & (mvirs < upperGroupMass))  # mask groups that pass the first mass cut
groupMasses = mvirs[groupMask] # to get all group masses in Msun instead of Msun/h
groupNumbers = subgr[inds[groupMask]]
firstSubhaloNumbers = inds[groupMask] 

#####################################
#####     Subhalo Mass Cuts     #####
#####################################
lowerSubhaloMass = 8 #8e10 Msun
upperSubhaloMass = 32 #3.2e11 Msun
lowerCurrentMass = 1 # 1e10 Msun

groupsList = list(firstSubhaloNumbers)
groups           = []
groupMasses      = []
groupRadii       = []
isoFlagHill      = []
pertGroupNum     = []
pertMassRatio    = []
pertSep          = []

# tidalIndex1p5Mpc = []  
# tidalIndex1Mpc   = []
# tidalIndex0p5Mpc = []
# isoFlag1p5Mpc    = []
# isoFlag1Mpc      = []
# isoFlag0p5Mpc    = []

subhaloMask = mvirs != 0 # to remove all groups with 0 mass
