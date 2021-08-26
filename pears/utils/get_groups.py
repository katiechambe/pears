""" 
Collects information about groups that pass group mass criteria

Dependencies:
-------------
ReadCats class from read_group_cats.py
SetupPaths class from paths.py

"""

__author__ = "Katie Chamberlain"
__status__ = "Beta -- forever."

import numpy as np
from utils.read_group_cats import ReadCats
from utils.paths import SetupPaths


class GetGroups:

    def __init__(
        self,
        snapshot=135,
        sim="Illustris",
        physics="dark",
        scale="dwarf",
        **kwargs
        ):

        """
        Defines bounds for group parameters and 

        Parameters:
        ----------_
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
        scale: str
            "dwarf" or "massive"
        kwargs: dict
            group_mass_min: float
                minimum mass of the group in 1e10Msun, default 8
            group_mass_max: float
                maximum mass of the group in 1e10Msun, default 50
            little_h: float
                definition of little h, default h=0.702
        """

        self.snapshot = snapshot
        self.sim = sim
        self.physics = physics
        self.scale = scale
        self.kwargs = kwargs
        
        self.group_min = self.kwargs.pop("group_mass_min", 8)
        self.group_max = self.kwargs.pop("group_mass_max", 50)
        self.h = self.kwargs.pop("little_h", 0.702)
        
        ReadCats.__init__(
            self, 
            snapshot=self.snapshot, 
            sim=self.sim, 
            physics=self.physics
            )

        self.masses = self.mvirs

        mask = self.masses > lowerGroupMass) & (self.masses < upperGroupMass)

        # mask = np.where((self.masses > lowerGroupMass) & (self.masses < upperGroupMass))




        print("got catalog")

        
        dict={{"group_mass_min":8},{"group_mass_max":50}}

        number_of_groups = len(self.mvirs)
        group_numbers = np.arange(0, number_of_groups,1)


groupMask =   # mask groups that pass the first mass cut
groupMasses = mvirs[groupMask] # to get all group masses in Msun instead of Msun/h
groupNumbers = subgr[inds[groupMask]]


subhaloMask = mvirs != 0 # to remove all groups with 0 mass

groupNum = subgr[i]





    def get_isolation_flag(self, group_number):
        xLo, xHi, yLo, yHi, zLo, zHi = subbox(groupPos[groupNum]) # in ckpc/h
        xGrps, yGrps, zGrps = groupPos[:,0], groupPos[:,1], groupPos[:,2]
        xPos, yPos, zPos = groupPos[groupNum][0], groupPos[groupNum][1], groupPos[groupNum][2]



    def save_groups(self):
        SetupPaths.__init__(self)
        save_path = self.path_groups + f"{self.sim}_{self.physics}_{self.scale}_{self.snapshot}.csv"

groupMass = mvirs[groupNum]
    groupRadius = rvirs[groupNum]
    groupMasses.append(groupMass) 
    groupRadii.append(groupRadius)   

    if xLo < xHi:
        xMask = np.logical_and(xGrps > xLo, xGrps < xHi)
    else:
        xMask = np.logical_or(xGrps > xLo, xGrps < xHi)

    if yLo < yHi:
        yMask = np.logical_and(yGrps > yLo, yGrps < yHi)
    else:
        yMask = np.logical_or(yGrps > yLo, yGrps < yHi)

    if zLo < zHi:
        zMask = np.logical_and(zGrps > zLo, zGrps < zHi)
    else: 
        zMask = np.logical_or(zGrps > zLo, zGrps < zHi)

    equivMask = np.logical_and(xGrps == xPos, yGrps == yPos) # to remove the group itself
    massMask  = (mvirs > groupMass) # only consider groups with masses higher than the group mass

    fullMask = xMask&yMask&zMask&~equivMask&subhaloMask&massMask

    if (~fullMask).all():
        massMask  = (mvirs > groupMass/10) # only consider groups with masses higher than the group mass
        fullMask = xMask&yMask&zMask&~equivMask&subhaloMask&massMask
        if (~fullMask).all():
            fullMask = xMask&yMask&zMask&~equivMask&subhaloMask



    linear = np.arange(0,len(fullMask),1)

    groupNumbers = linear[fullMask]
    masses       = mvirs[fullMask]
    massRatios    = masses/groupMass
    radiuses     = rvirs[fullMask]

    xPoses = groupPos[fullMask][:,0]
    yPoses = groupPos[fullMask][:,1]
    zPoses = groupPos[fullMask][:,2]

    primaryPos = groupPos[fullMask]
    numberGroups = len(primaryPos)
    secondaryPos = np.column_stack([xPos*np.ones(numberGroups),yPos*np.ones(numberGroups),zPos*np.ones(numberGroups)])

    correctedSeparations = np.array([np.linalg.norm(i) for i in np.array(vector(primaryPos,secondaryPos,scale))]) # physical kpc

    hillRadii = hillRadius(groupMass, correctedSeparations, masses)

    if (hillRadii < groupRadius).any():
        isolation_flag = False
    else:
        isolation_flag = True

        return isolation_flag 



    
    pertInd = np.where(np.min(hillRadius) == hillRadius)[0][0]
    perturber = groupNumbers[pertInd]
    pertGroupNum.append(perturber)
    pertMassRatio.append( (mvirs[perturber]) / groupMass )
    pertSep.append( correctedSeparations[pertInd] )

    count = groupsList.index(i)
    if count%100 == 0:
        print('Have gone through ',count," groups out of ", len(firstSubhaloNumbers)+1)
print('got groups')

#######################################
# collecting the groups and exporting #
#######################################
# groupNumbers = np.array(groups)
zipped = list(zip(groups, groupMasses, groupRadii, isoFlagHill, pertGroupNum, pertMassRatio, pertSep))
df = pd.DataFrame(data = zipped, columns=['Group Numbers','Group Mass','Group Radius','Isolated Flag Hill','Perturbing Group Number','Perturbing Group Mass Ratio','Perturbing Group Separation'])
# zipped = list(zip(groups, groupMasses, isoFlag0p5Mpc, tidalIndex0p5Mpc, isoFlag1Mpc, tidalIndex1Mpc, isoFlag1p5Mpc, tidalIndex1p5Mpc, pertGroupNum, pertMassRatio, pertSep))
# df = pd.DataFrame(data = zipped, columns=['Group Numbers','Group Mass','Isolated 0.5Mpc Flag','Tidal Index 0.5Mpc','Isolated 1Mpc Flag','Tidal Index 1Mpc', 'Isolated 1.5Mpc Flag','Tidal Index 1.5Mpc','Perturbing Group Number','Perturbing Group Mass Ratio','Perturbing Group Separation'])
# note, group mass is in 10^10 Msun
df.to_csv(saveFile,index=False,header=True)
