""" 
Dependencies:
-------------

"""
__author__ = "Katie Chamberlain"
__status__ = "Beta - forever~"
__date__   = "July 2022"

import sys
import h5py
import numpy as np
from utils.paths import SetupPaths

paths = SetupPaths()
  
    
class compile_summary:
    
    def __init__(
        self,
        sim="Illustris", 
        phys="dark",
        size="dwarf",
        subset="major",
        **kwargs
        ):
    
        """
        Collects the 

        Parameters:
        -----------
        sim: str
            "Illustris" or "TNG"
            to specify which simulation
        physics: str
            "dark" or "hydro"
            to specify which simulation
        size: str
            "dwarf" or "massive"
        subset: str
            to specify the subset
            choose from: major, minor, all
        """
        
        self.sim = sim
        self.phys = phys
        self.size = size
        self.subset = subset
        self.kwargs = kwargs
        self.redcutoff = self.kwargs.pop("redshift_cutoff", 1000)
        

    def get_data(self):
        
        summarydata = {"Redshift":[],
                            "Snapshot":[],
                            "Number pairs":[],
                            "Number primaries":[],
                           # "Number tertiaries":[],
                            "Ratio pairs":[],
                            "Median Separation": [],
                            "Median Separation Quartiles": [],
                            "Mean Separation": [],
                            "Mean Separation Std": [],
                            "Median RelVel": [],
                            "Median RelVel Quartiles": [],
                            "Mean RelVel": [],
                            "Mean RelVel Std": [],
                            "Mean Lowsep Counts": [],
                            "Mean Lowsep Counts Std": []}
#         summarydata = {}
        
        self.snapshots = {"Illustris":np.arange(0,136), 
                     "TNG":np.arange(0,100)}
        
        for snapshot in self.snapshots[self.sim]:
            pair_path = f"{self.sim}_{snapshot}_10.hdf5"
            pair_data = h5py.File(f"{paths.path_pairs}{pair_path}", "r")
      

            try:
                data = pair_data[self.phys]
                red_coor = pair_data['Header'].attrs['Redshift']

            except KeyError:
                print(f"{snapshot} does not exist -- skipping")

            else:
                if red_coor < self.redcutoff:
                    med_mask = np.array(data["Realization"]) == -1

                    primary_stell = np.array(data["Sub1 Stellar Mass"])

                    # primary stellar mass will be different for dwarf and massive pairs!~
                    if self.size == "dwarf":
                        primary_analog = (primary_stell > 0.01) & (primary_stell < 0.5)
                    elif self.size == "massive":
                        primary_analog = (primary_stell > 0.5) & (primary_stell < 10) # this doesn't matter, can be changed! 

                    # get major pairs only  
                    majors = np.array(data["Stellar Mass Ratio"]) > 1/4
                    minors = (np.array(data["Stellar Mass Ratio"]) < 1/4) & (np.array(data["Stellar Mass Ratio"]) > 1/10)
                    pair = majors+minors
                    
                    minsep = np.array(data["Separation"]) >10
                    lowsep = np.array(data["Separation"]) < 50
                    highsep = np.array(data["Separation"]) > 50
                    
                    lowvel = np.array(data["RelVel"]) < 100
                    highvel = np.array(data["RelVel"]) > 100
                    
                    countarray = []
                    
                    for i in np.arange(-1,10):
                        realmask = np.array(data["Realization"]) == i
                        counts = np.sum(lowsep & primary_analog & pair & minsep & realmask)
                        numprims = np.sum(primary_analog & realmask)
                        countarray.append(counts/numprims)

                    masks = {"primaries":primary_analog,
                             "allpairs":primary_analog & minsep,  # this should be & pairs
                             "median": primary_analog & med_mask,
                             "major":majors & primary_analog & minsep, 
                             "minor":minors & primary_analog & minsep,
                             "majormed":majors & med_mask & primary_analog & minsep,
                             "minormed":minors & med_mask & primary_analog & minsep,
                             "lowsep":lowsep & primary_analog & pair & minsep,
                             "highsep":highsep & primary_analog & pair & minsep,
                             "lowsepmajor":lowsep & primary_analog & pair & minsep & majors,
                             "highsepmajor":highsep & primary_analog & pair & minsep & majors,
                             "lowvel":lowvel & primary_analog & pair & minsep,
                             "highvel":highvel & primary_analog & pair & minsep,
                             "lowvelmajor":lowvel & primary_analog & pair & minsep & majors,
                             "highvelmajor":highvel & primary_analog & pair & minsep & majors}

                    use_mask = masks[self.subset]

                    if all(use_mask == False):
                        print(f"{snapshot} {self.phys} {self.size} no pairs in subset {self.subset}")
                        continue

                    sublist = {} # dictionary  of subhalos in group
                    for key, val in data.items():
                        sublist[key] = np.array(val)[use_mask]

                    # redshift and count info
                    redshift = red_coor
                    num_pairs = np.sum(use_mask)
                    
                    # put tertiaries here maybe?
                    num_primaries = np.sum(primary_analog)
                    ratio_pairs = num_pairs/num_primaries

                    # separation info
                    med_sep = np.median(sublist['Separation'])
                    qs_sep = np.percentile(sublist['Separation'],[16,84])
                    mean_sep = np.mean(sublist['Separation'])
                    std_sep = np.std(sublist['Separation'])

                    # velocity info
                    med_vel = np.median(sublist['RelVel'])
                    qs_vel = np.percentile(sublist['RelVel'],[16,84])
                    mean_vel = np.mean(sublist['RelVel'])
                    std_vel = np.std(sublist['RelVel'])
                    
                    med_counts = np.mean(countarray)
                    std_counts = np.std(countarray)
                    
                    summarylist = {"Redshift":redshift,
                            "Snapshot":snapshot,
                            "Number pairs":num_pairs,
                            "Number primaries":num_primaries,
                            "Ratio pairs":ratio_pairs,
                            "Median Separation": med_sep,
                            "Median Separation Quartiles": qs_sep,
                            "Mean Separation": mean_sep,
                            "Mean Separation Std": std_sep,
                            "Median RelVel": med_vel,
                            "Median RelVel Quartiles": qs_vel,
                            "Mean RelVel": mean_vel,
                            "Mean RelVel Std": std_vel,
                            "Mean Lowsep Counts": med_counts,
                            "Mean Lowsep Counts Std": std_counts}
                    
#                     summarylist = np.array([redshift, snapshot, 
#                                          num_pairs, num_primaries, ratio_pairs, 
#                                          med_sep, qs_sep, mean_sep, std_sep,
#                                          med_vel, qs_vel, mean_vel, std_vel,
#                                          med_counts, std_counts])

                    for ind, key in enumerate(summarylist.keys()):
#                         if snapshot == self.snapshots[self.sim][0]:
#                             summarydata[key]=summarylist[key]
#                         else:
                        summarydata[key].append(summarylist[key])

        return summarydata
    
    def get_snapshot(self, snapnum):
        
        pair_path = f"{self.sim}_{snapnum}_10.hdf5"
        pair_data = h5py.File(f"{paths.path_pairs}{pair_path}", "r")

        try:
            data = pair_data[self.phys]
            red_coor = pair_data['Header'].attrs['Redshift']

        except KeyError:
            print(f"snapshot {snapshot}: data does not exist - try another snapshot")

        else:
            med_mask = np.array(data["Realization"]) == -1

            primary_stell = np.array(data["Sub1 Stellar Mass"])

            # primary stellar mass will be different for dwarf and massive pairs!~
            if self.size == "dwarf":
                primary_analog = (primary_stell > 0.01) & (primary_stell < 0.5)
            elif self.size == "massive":
                primary_analog = (primary_stell > 0.5) & (primary_stell < 10) # this doesn't matter, can be changed! 

            # get major pairs only  
            majors = np.array(data["Stellar Mass Ratio"]) > 1/4
            minors = (np.array(data["Stellar Mass Ratio"]) < 1/4) & (np.array(data["Stellar Mass Ratio"]) > 1/10)
            pair = majors+minors

            minsep = np.array(data["Separation"]) >10
            lowsep = np.array(data["Separation"]) < 50
            highsep = np.array(data["Separation"]) > 50
            
            lowvel = np.array(data["RelVel"]) < 100
            highvel = np.array(data["RelVel"]) > 100

            countarray = []

            for i in np.arange(-1,10):
                realmask = np.array(data["Realization"]) == i
                counts = np.sum(lowsep & primary_analog & pair & minsep & realmask)
                numprims = np.sum(primary_analog & realmask)
                countarray.append(counts/numprims)

            masks = {"primaries":primary_analog,
                     "allpairs":primary_analog & minsep,
                     "median": primary_analog & med_mask,
                     "major":majors & primary_analog & minsep, 
                     "minor":minors & primary_analog & minsep,
                     "majormed":majors & med_mask & primary_analog & minsep,
                     "minormed":minors & med_mask & primary_analog & minsep,
                     "lowsep":lowsep & primary_analog & pair & minsep,
                     "highsep":highsep & primary_analog & pair & minsep,
                     "lowsepmajor":lowsep & primary_analog & pair & minsep & majors,
                     "highsepmajor":highsep & primary_analog & pair & minsep & majors,
                     "lowvel":lowvel & primary_analog & pair & minsep,
                     "highvel":highvel & primary_analog & pair & minsep,
                     "lowvelmajor":lowvel & primary_analog & pair & minsep & majors,
                     "highvelmajor":highvel & primary_analog & pair & minsep & majors}

            use_mask = masks[self.subset]
            num_pairs = np.sum(use_mask)
            num_primaries = np.sum(primary_analog)
            ratio_pairs = num_pairs/num_primaries

            if all(use_mask == False):
                print(f"{snapshot} {self.phys} {self.size} no pairs in subset {self.subset}")
                           
            # dictionary  of subhalos in group
            sublist = {"Snapshot":snapnum,
                       "Redshift":red_coor,
                       "Simulation":self.sim,
                       "Number Pairs":np.sum(use_mask),
                       "Number Primaries":np.sum(primary_analog),
                       "Number Tertiaries":np.sum(primary_analog),
                       "Ratio Pairs": num_pairs/num_primaries}
            
            for key, val in data.items():
                sublist[key] = np.array(val)[use_mask]

        return sublist
    
