""" 
Collects information about groups that pass group mass criteria

Dependencies:
-------------
ReadCats class from read_group_cats.py

"""
__author__ = "Katie Chamberlain"
__status__ = "Beta - forever~"
__date__   = "September 2021"

import numpy as np
from utils.read_group_cats import ReadCats
from astropy.table import QTable
import astropy.units as u

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
            omega_m:
            omega_r:
            omega_l:

        """

        self.snapshot = snapshot
        self.sim = sim
        self.physics = physics
        self.size = size
        self.kwargs = kwargs
        
        # set mass scale for dwarf and massive groups
        if self.size == "dwarf":
            self.group_min = 8
            self.group_max = 50

        elif self.size == "massive":
            self.group_min = 100
            self.group_max = 400

        self.group_min = self.kwargs.pop("group_mass_min", self.group_min)
        self.group_max = self.kwargs.pop("group_mass_max", self.group_max)
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
        self.submass_phys = self.submass/self.h 
        self.subpos_phys = self.subpos*self.scale/self.h 
        low_cut = (self.mvirs_phys > self.group_min)
        high_cut = (self.mvirs_phys < self.group_max)
        non_zero = self.mvirs != 0 
        mass_cut = low_cut&high_cut&non_zero

        # find the Group ID #s of all groups that pass the cut
        self.pass_numbers = np.where(mass_cut==True)[0]
        self.pass_mvir = self.mvirs_phys[self.pass_numbers]
        self.pass_rvir = self.rvirs_phys[self.pass_numbers]
        self.pass_nsubs = self.nsubs[self.pass_numbers]

        # find intersection of all halos in groups that pass mass cut
        # and that have current mass >1e9Msun
        subhalo_group_mask = np.in1d(self.subgr, self.pass_numbers)
        subhalo_mass_mask = self.submass_phys >= 0.1 
        self.subhalo_mask = np.where(subhalo_group_mask&subhalo_mass_mask)
    
    def save_groups(self):
        '''
        Save the groups that pass the mass cuts
        '''
        save_path = f"{self.sim}_{self.physics}_{self.size}_{self.snapshot}.ecsv"

        group_number = u.Quantity(self.pass_numbers,
                                  dtype=np.int64)
        group_number.info.description = "Group Number"

        group_mass = self.pass_mvir * u.Unit(1e10 * u.Msun)
        group_mass.info.description = "Physical mass from Group_M_TopHat200"

        group_radius = self.pass_rvir * u.kpc
        group_radius.info.description = "Physical radius from Group_R_TopHat200"

        group_nsubs = u.Quantity(self.pass_nsubs,
                                 dtype=np.int64)
        group_nsubs.info.description = "Number of subhalos in group"

        t = QTable()
        t['group_number'] = group_number
        t['group_mass'] = group_mass
        t['group_radius'] = group_radius
        t['group_nsubs'] = group_nsubs
        t.write(self.path_groups + save_path,
                overwrite=True)

        print(f"Saved groups at {self.sim}_{self.physics}_{self.size}_{self.snapshot}.csv")

# TODO: fix this up! 
# TODO: make sure to save the redshift of each snapshoit in the metadata of the ecsv file! 
    def save_subhalos(self):
            '''
            Save the subhalo data
            '''
                # max mass and corresponding snapnumber
            inst = mergerClass.MergerTree(snap, subID)
                maxMass, maxMassSnap    = inst.maxMass()[0]/h, inst.maxMass()[1]
                # i is group number, subID is subhalo ID, 
                data.append([groupID, subID, massz0, maxMass, maxMassSnap, ])
            


            save_path = f"{self.sim}_{self.physics}_{self.size}_{self.snapshot}.ecsv"

            subhalo_number = u.Quantity(self.pass_numbers,
                                    dtype=np.int64)
            subhalo_number.info.description = "Group Number"

            subhalo_mass = self.pass_mvir * u.Unit(1e10 * u.Msun)
            subhalo_mass.info.description = "Physical mass from Group_M_TopHat200"

            subhalo_nsubs = u.Quantity(self.pass_nsubs,
                                    dtype=np.int64)
            subhalo_nsubs.info.description = "Number of subhalos in group"

            t.meta?? 


    subhalo_pos, subhalo_vel, subhalo_mass, subhalo_grnmb, subhalo_IDs, subhalo_maxmass, subhalo_maxmass_snap


            t = QTable()
            t['group_number'] = group_number
            t['group_mass'] = group_mass
            t['group_radius'] = group_radius
            t['group_nsubs'] = group_nsubs
            t.write(self.path_subhalos + save_path,
                    overwrite=True)

            print(f"Saved groups at {self.sim}_{self.physics}_{self.size}_{self.snapshot}.csv")
