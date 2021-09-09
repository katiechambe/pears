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
from astropy.table import QTable
import astropy.units as u
# import pandas as pd


class GetGroups:

    def __init__(
        self,
        snapshot=135,
        sim="Illustris",
        physics="dark",
        size="dwarf",
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
            "Illustris" or "TNG"
            to specify which simulation
        physics: str
            "dark" or "hydro"
            to specify which simulation
        size: str
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
        self.size = size
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
        self.scale = 1/(1+self.redshift)

        # to pull out all of the groups within the bounds
        # Note: virial mass and virial radius are both in physical units
        self.mvirs_phys = self.mvirs/self.h 
        self.rvirs_phys = self.rvirs*self.scale/self.h
        low_cut = (self.mvirs_phys > self.group_min)
        high_cut = (self.mvirs_phys < self.group_max)
        non_zero = self.mvirs != 0 
        mass_cut =  low_cut&high_cut&non_zero

        # find the Group ID #s of all groups that pass the cut
        self.pass_numbers = np.where(mass_cut==True)[0]
        self.pass_mvir = self.mvirs_phys[self.pass_numbers]
        self.pass_rvir = self.rvirs_phys[self.pass_numbers]
        self.pass_nsubs = self.nsubs[self.pass_numbers]

    def save_groups(self):
        SetupPaths.__init__(self)
        save_path = f"{self.sim}_{self.physics}_{self.size}_{self.snapshot}.ecsv"

        group_number = u.Quantity(self.pass_numbers,dtype=np.int64)
        group_number.info.description = "Group Number"

        group_mass = self.pass_mvir*u.Unit(1e10 * u.Msun)
        group_mass.info.description = "Physical mass from Group_M_TopHat200"

        group_radius = self.pass_rvir*u.kpc
        group_radius.info.description = "Physical radius from Group_R_TopHat200"

        group_nsubs = u.Quantity(self.pass_nsubs,dtype=np.int64)
        group_nsubs.info.description = "Number of subhalos in group"

        t = QTable()
        t['group_number'] = group_number
        t['group_mass'] = group_mass
        t['group_radius'] = 
        t['group_nsubs'] = 
        self.pass_numbers.info.description
        t.write(self.path_groups + save_path)


        # zipped = list(zip(self.pass_numbers, 
        #                   self.pass_mvir*u.Unit(1e10 * u.Msun), 
        #                   self.pass_rvir*u.kpc, 
        #                   self.pass_nsubs))

        # df = pd.DataFrame(data=zipped, 
        #     columns=["Group Number","Group Mass","",""])

        df = pd.DataFrame(data = )

        print(f"Saved groups at {self.sim}_{self.physics}_{self.size}_{self.snapshot}.csv")


print('got groups')

#######################################
# collecting the groups and exporting #
#######################################
# groupNumbers = np.array(groups)
df = pd.DataFrame(data = zipped, columns=['Group Numbers','Group Mass','Group Radius','Isolated Flag Hill','Perturbing Group Number','Perturbing Group Mass Ratio','Perturbing Group Separation'])
# zipped = list(zip(groups, groupMasses, isoFlag0p5Mpc, tidalIndex0p5Mpc, isoFlag1Mpc, tidalIndex1Mpc, isoFlag1p5Mpc, tidalIndex1p5Mpc, pertGroupNum, pertMassRatio, pertSep))
# df = pd.DataFrame(data = zipped, columns=['Group Numbers','Group Mass','Isolated 0.5Mpc Flag','Tidal Index 0.5Mpc','Isolated 1Mpc Flag','Tidal Index 1Mpc', 'Isolated 1.5Mpc Flag','Tidal Index 1.5Mpc','Perturbing Group Number','Perturbing Group Mass Ratio','Perturbing Group Separation'])
# note, group mass is in 10^10 Msun
df.to_csv(saveFile,index=False,header=True)
