import h5py
import os.path
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gd
from utils.paths import SetupPaths
import matplotlib.patches as  mpatches
from matplotlib.ticker import FormatStrFormatter
from utils.get_summary_data import compile_summary
from utils.vectorCorrection import vectorCorrection as vector


paths = SetupPaths()

# defining equivalent redshifts
zs = {"z":np.array([0,1,2,3,4]), 
#       "zill":np.array([135,85,68,60,56]), 
      "ztng":np.array([99,50,33,25,21])}

# make functions to get data at the requested snapshot
def get_primmask(primstells, size):
    if size == "dwarf":
        mask = (primstells > 0.01) & (primstells < 0.5)
    elif size == "massive":
        mask = (primstells > 0.5) & (primstells < 10)
    return mask

def get_groupmask(groupmass, size):
    if size == "dwarf":
        mask = (groupmass > 8) & (groupmass < 50)
    elif size == "massive":
        mask = (groupmass > 100) & (groupmass < 650)
    return mask

class EmptyFile(Exception): pass
class SkipRedshift(Exception): pass


def get_sepvel(reals, key, scaled=False, errorprint=False, redshiftcutoff=True): 
    snapshots = np.arange(0,100,1)
    snapshots = np.delete(snapshots, np.where(snapshots==48)[0])
    redcutoff = 4.2
        
    redshifts = []  
    med_majdw, med_mindw, med_majma, med_minma = [], [], [], []
    med_majdiff, med_mindiff, med_ddiff, med_mdiff = [], [], [], []
    quart_majdw, quart_mindw, quart_majma, quart_minma = [], [], [], []
    quart_majdiff, quart_mindiff, quart_ddiff, quart_mdiff  = [], [], [], []
            
    for snap in snapshots:  
        try:
            pair_path = f"TNG_{snap}_{reals}.hdf5"
            pair_data = h5py.File(f"{paths.path_pairs}{pair_path}", "r")
            
            if np.size(pair_data) == 0:
                raise EmptyFile
                
            redshift = pair_data['Header'].attrs['Redshift']
            
            if redshiftcutoff & ( redshift > redcutoff) :
                raise SkipRedshift
                
            if (len(pair_data['pairs']["hydro"]['Group ID']) == 0):    
                raise EmptyFile
                
            pair = pair_data["pairs"]["hydro"]
            priStell = np.array(pair["Sub1 Stellar Mass"])
            secStell = np.array(pair["Sub2 Stellar Mass"])
            pairGroups = np.array(pair["Group Mass"])
            pairGrRads = np.array(pair["Group Radius"])
            pairReals = np.array(pair["Realization"])
            seps = np.array(pair["Separation"]) 
            vels = np.array(pair["RelVel"])
                        
            majors = (secStell/priStell > 1/4)
            minors = (secStell/priStell > 1/10) & (secStell/priStell < 1/4)
            pair_lowsep = (seps > 10)
            
            if key == "Separation":
                scaleddat = seps / pairGrRads
            elif key == "RelVel":
                G = 4.3009173e4 # in km^2 kpc / (1e10M⊙ s^2)
                vvir = np.sqrt(G*pairGroups / pairGrRads)
                scaleddat = vels / (vvir)  
                
            if scaled: 
                dat = scaleddat
            elif key == "Separation":
                dat = seps
            elif key == "RelVel":
                dat = vels
            else:
                dat = np.array(pairs[key]) 
            
                 ## dwarfs
            pair_pri_dwarf = get_primmask(priStell, "dwarf")
            pair_group_dwarf = get_groupmask(pairGroups, "dwarf")
            
                # defining combined masks 
            pair_mask_dwarf = pair_pri_dwarf & pair_group_dwarf & pair_lowsep
            
            ## massive
            pair_pri_massive = get_primmask(priStell, "massive")
            pair_group_massive = get_groupmask(pairGroups, "massive")

                # defining combined masks 
            pair_mask_massive = pair_pri_massive & pair_group_massive & pair_lowsep
            
            real_majdw = []
            real_mindw = []
            real_majma = []
            real_minma = []
            real_majdiff = []
            real_mindiff = []
            real_ddiff = []
            real_mdiff = []
            
            realizations = np.unique( pairReals )

            for real in realizations:
                pair_real = pairReals == real

                mask_majdw = pair_real & pair_mask_dwarf & majors
                mask_mindw = pair_real & pair_mask_dwarf & minors
                mask_majma = pair_real & pair_mask_massive & majors
                mask_minma = pair_real & pair_mask_massive & minors

                majdw_xx = np.median( dat[mask_majdw] )
                mindw_xx = np.median( dat[mask_mindw] )
                majma_xx = np.median( dat[mask_majma] )
                minma_xx = np.median( dat[mask_minma] )

                real_majdw.append( majdw_xx )
                real_mindw.append( mindw_xx )
                real_majma.append( majma_xx )
                real_minma.append( minma_xx )
                real_majdiff.append( majma_xx - majdw_xx)
                real_mindiff.append( minma_xx - mindw_xx)
                real_ddiff.append( majdw_xx - mindw_xx)
                real_mdiff.append( majma_xx - minma_xx)

            lower, upper = 16,84         
            redshifts.append( redshift )
            
            med_majdw.append( np.median(real_majdw) )
            med_mindw.append( np.median(real_mindw) )
            med_majma.append( np.median(real_majma) )
            med_minma.append( np.median(real_minma) )
            med_majdiff.append( np.median(real_majdiff) )
            med_mindiff.append( np.median(real_mindiff) )
            med_ddiff.append( np.median(real_ddiff) )
            med_mdiff.append( np.median(real_mdiff) )
            
            quart_majdw.append( np.percentile( real_majdw, [lower,upper] ) )
            quart_mindw.append( np.percentile( real_mindw, [lower,upper] ) )
            quart_majma.append( np.percentile( real_majma, [lower,upper] ) )
            quart_minma.append( np.percentile( real_minma, [lower,upper] ) )
            quart_majdiff.append( np.percentile( real_majdiff, [lower,upper] ) )
            quart_mindiff.append( np.percentile( real_mindiff, [lower,upper] ) )
            quart_ddiff.append( np.percentile( real_ddiff, [lower,upper] ) )
            quart_mdiff.append( np.percentile( real_mdiff, [lower,upper] ) )

        except KeyError:
            if errorprint: print(f'skipping {snap} for KeyError. Please debug')
            continue
            
        except EmptyFile:
            if errorprint: print(f"skipping {snap}, empty file")
            continue
            
        except SkipRedshift:
            if errorprint: print(f"skipping {snap}, redshift out of range")
                
    scaled_dictionary = {"z": np.array(redshifts),

                        "Median Major Dwarf": np.array(med_majdw),
                        "Median Minor Dwarf": np.array(med_mindw),
                        "Median Major Massive": np.array(med_majma),
                        "Median Minor Massive": np.array(med_minma),
                        "Median Major Difference": np.array(med_majdiff),
                        "Median Minor Difference": np.array(med_mindiff),
                        "Median Dwarf Difference": np.array(med_ddiff),
                        "Median Massive Difference": np.array(med_mdiff),

                        "Quartile Major Dwarf": np.array(quart_majdw),
                        "Quartile Minor Dwarf": np.array(quart_mindw),
                        "Quartile Major Massive": np.array(quart_majma),
                        "Quartile Minor Massive": np.array(quart_minma),
                        "Quartile Major Difference": np.array(quart_majdiff),
                        "Quartile Minor Difference": np.array(quart_mindiff),
                        "Quartile Dwarf Difference": np.array(quart_ddiff),
                        "Quartile Massive Difference": np.array(quart_mdiff)}
    
    return scaled_dictionary

def make_sepveldata(reals):
        #check if files exists
    filepath = f"{paths.path_plotdata}sepvel.hdf5"
    if not os.path.isfile(filepath):
        print("file does not exist...")
        print("creating file")
        f = h5py.File(filepath, 'w')
        print("file created successfully. adding header...")
        header_dict = {"1000 Reals - Quartile Range":"16-84%",
                       "Simulation":"TNG100-1 (Hydro)"}

        dset = f.create_group('/Header')
        for key in header_dict.keys():
            dset.attrs[key] = header_dict[key]
            
        print("header added successfully")
    else:
        print("file exists...")
        f = h5py.File(filepath, 'r+')
        
        
    print("checking to see if data exists for this number of realizations")
    
    if f.get(f"{reals} Realizations") is not None:
        print("data already exists!")
        f.close()
              
    else:
        print("data does not exist...")
        print("creating data tables...")
        
        pairs = [("Separation",True),("Separation",False),("RelVel", True),("RelVel", False)]
        labels = ["Scaled Separation", "Separation", "Scaled Velocity", "Velocity"]         
        
        for pa, la in zip(pairs, labels):
            datums = get_sepvel(reals, pa[0], scaled=pa[1])
            print(f"collected data for {pa}")

            for key, val in datums.items():
                val = np.array(val)
                dset = f.create_dataset(f'/{reals} Realizations/{la}/{key}', 
                                        shape=val.shape,
                                        dtype=val.dtype)
                dset[:] = val
                
        f.close()
        print("data saved~")
    
    
    
print("running make sepvel data")
make_sepveldata(10)
print("done with 10 reals")
make_sepveldata(100)
print("done with 100 reals")
make_sepveldata(1000)
print("finally done with 1000 reals")

