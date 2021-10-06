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
from utils.merger_trees import TraceMergerTree
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
        -----------
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

        if self.sim == "Illustris":
            self.little_h = 0.704

        elif self.sim == "TNG":
            self.little_h = 0.6774
        
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
        self.subhalo_ids = np.where(subhalo_group_mask&subhalo_mass_mask)[0]
        self.subhalo_masses = self.submass_phys[self.subhalo_ids]

        self.save_path = f"{self.sim}_{self.physics}_{self.size}_{self.snapshot}.ecsv"
    


# TODO: fix this up! 
    def save_subhalos(self):
        '''
        Save the subhalo data for groups that pass mass cuts.
        Note:
        -----
        This may be slow because it has to look up the merger trees for 
        each of the subhalos 
        '''
        x = []
        mask = []
        prev_group=None
        for subid in self.subhalo_ids: 
            current_group = self.subgr[subid]



# dataa = h5py.File(f"{groups.path_data}max_masses/{sim}_{phys}.hdf5", 'r')


            # only pay attention to first 5 subhalos
            if prev_group is None:
                prev_group = current_group     
                i=-1
            elif prev_group==current_group:
                i+=1
            elif prev_group!=current_group:
                i=0

            if i < 5:
                try:
                    x.append(TraceMergerTree(
                            snapshot=hey.snapshot,
                            subfindID=subid,
                            sim=self.sim, 
                            physics=self.physics).maxmass )
                    mask.append(True)

                except AttributeError:
                    print(f'Could not find merger tree for {sim} {physics} {subid}')
                    mask.append(False)

            else:
                mask.append(False)
            prev_group=current_group

        x = np.array(x)

        sub_max = x[:,0]
        sub_max_snap = x[:,1]

        subhalo_groupnum = u.Quantity( self.subgr[self.subhalo_ids][mask], dtype=np.int64)
        subhalo_groupnum.info.description = "Group Number"

        subhalo_id = u.Quantity( self.subhalo_ids[mask], dtype=np.int64) 
        subhalo_id.info.description = "Subhalo ID"

        subhalo_mass =  self.submass_phys[self.subhalo_ids][mask] * u.Unit(1e10 * u.Msun)
        subhalo_mass.info.description = "Physical mass bound to subhalo"

        subhalo_pos = self.subpos_phys[self.subhalo_ids][mask] * u.kpc
        subhalo_pos.info.description = "Physical position of subhalo in the box"

        subhalo_vel = self.subvel[self.subhalo_ids][mask] * u.km / u.s
        subhalo_vel.info.description = "Peculiar velocity of subhalo"

        subhalo_maxmass =  sub_max * u.Unit(1e10 * u.Msun)
        subhalo_maxmass.info.description = "Maximum mass (physical) ever achieved"

        subhalo_maxmass_snap = u.Quantity( sub_max_snap , dtype=np.int64)
        subhalo_maxmass_snap.info.description = "Snapshot where maximum mass is achieved"
            
        t = QTable(meta={"snapshot":self.snapshot, "redshift":self.redshift})

        t['subhalo_groupnum'] = subhalo_groupnum
        t['subhalo_id'] = subhalo_id
        t['subhalo_mass'] = subhalo_mass
        t['subhalo_pos'] = subhalo_pos
        t['subhalo_vel'] = subhalo_vel
        t['subhalo_maxmass'] = subhalo_maxmass
        t['subhalo_maxmass_snap'] = subhalo_maxmass_snap

        t.write(self.path_subhalos + self.save_path,
                overwrite=True)

        print(f"Saved subhalos at {self.sim}_{self.physics}_{self.size}_{self.snapshot}.ecsv")


    def save_groups(self):
        '''
        Save the groups that pass the mass cuts
        '''
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

        t = QTable(meta={"snapshot":self.snapshot, "redshift":self.redshift})
        t['group_number'] = group_number
        t['group_mass'] = group_mass
        t['group_radius'] = group_radius
        t['group_nsubs'] = group_nsubs
        t.write(self.path_groups + self.save_path,
                overwrite=True)

        print(f"Saved groups at data/groups/{self.sim}_{self.physics}_{self.size}_{self.snapshot}.ecsv")
