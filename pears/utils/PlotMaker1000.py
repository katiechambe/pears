'''
 PlotMaker1000 makes your favorite plot! 
 Just pick a plot name below ~                   
'''

Plot2Make = "counts"
Plot2Make = "pairfrac"
# Plot2Make = "sep"
# Plot2Make = "vel"
# Plot2Make = "sepdist"
# Plot2Make = "veldist"
# Plot2Make = "scaledsep"
# Plot2Make = "scaledvel"

maketype = "plot"
maketype = "data"
maketype = "both"

numreals = 1000

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
__author__ = "Katie Chamberlain"
__status__ = "Beta - forever~"
__date__   = "Feb 2023 - unedited "

import os.path
import h5py
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.gridspec as gd
from utils.paths import SetupPaths
import matplotlib.patches as  mpatches
from matplotlib.ticker import FormatStrFormatter
from utils.get_summary_data import compile_summary
from utils.vectorCorrection import vectorCorrection as vector


paths = SetupPaths()

plt.show();
plt.rcParams.update({'font.size':20,"xtick.direction":"in","ytick.direction":"in", 
                     "xtick.top":True, "ytick.right":True,"text.usetex":False,
                     "xtick.labelsize":18,"ytick.labelsize":18})

# defining color palette for plotting
palette = {"Illustris dark": "#009292", "Illustris hydro": "#B6DAFF",
           "TNG dark": "#930200", "TNG hydro": "#FFB5DC",
           "dwarf":"olive","massive":"salmon", "difference":"#2C1D11", "difference2":"#464646"}

alphas = {"maj": 0.7, "min": 0.3}

# defining equivalent redshifts
zs = {"z":np.array([0,1,2,3,4]), 
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

def get_counts( size, reals, errorprint=False, redshiftcutoff=True):    
    snapshots = np.arange(0,100,1)
    snapshots = np.delete(snapshots, np.where(snapshots==48)[0])
    redcutoff = 4.2
        
    redshifts = []
    medone, medtwo, medtot, medmaj, medmin, medpair, medmajfrac, medminfrac, medtotfrac = [], [], [], [], [], [], [], [], []
    quartone, quarttwo, quarttot, quartmaj, quartmin, quartpair, quartmajfrac, quartminfrac, quarttotfrac = [], [], [], [], [], [], [], [], []

    for snap in snapshots:  
        singleprims, doubleprims, totalprims = [], [], []
        majorpairs, minorpairs, totalpairs = [], [], []
        majfrac, minfrac, totfrac = [], [], []
        
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
                
            unpair = pair_data["unpaired"]["hydro"]
            unpairStells = np.array(unpair["Sub1 Stellar Mass"])
            unpairGroups = np.array(unpair["Group Mass"])
            unpairReals = np.array(unpair['Realization'])
            
            pair = pair_data["pairs"]["hydro"]
            priStell = np.array(pair["Sub1 Stellar Mass"])
            secStell = np.array(pair["Sub2 Stellar Mass"])
            pairGroups = np.array(pair["Group Mass"])
            pairReals = np.array(pair["Realization"])
            seps = np.array(pair["Separation"]) 
            
            # subset masks for unpaired
            unpair_pri = get_primmask(unpairStells, size)
            unpair_group = get_groupmask(unpairGroups, size)                
            
            pair_pri = get_primmask(priStell, size)
            pair_group = get_groupmask(pairGroups, size)

            majors = (secStell/priStell > 1/4)
            minors = (secStell/priStell > 1/10) & (secStell/priStell < 1/4)
            allpairs = (majors + minors)
            pair_lowsep = (seps > 10)

            # defining combined masks 
            unpair_mask = unpair_pri & unpair_group
            primary_mask = pair_pri & pair_group
            pair_mask = pair_pri & pair_group & pair_lowsep & allpairs
                                                           
            for real in np.unique(unpairReals):                  
                # make realization masks
                unpair_real = unpairReals == real
                pair_real = pairReals == real

                # make count values for single realization
                numone = np.count_nonzero(unpair_mask & unpair_real)
                numtwo = np.count_nonzero(primary_mask & pair_real)
                numtot = numone + numtwo
                nummaj = np.count_nonzero(pair_mask & pair_real & majors)
                nummin = np.count_nonzero(pair_mask & pair_real & minors)
                numpair = np.count_nonzero(pair_mask & pair_real)
                
                if numtot == 0:
                    continue
                
                # collect count vals for all reals
                singleprims.append(numone)
                doubleprims.append(numtwo)
                totalprims.append(numtot)
                majorpairs.append(nummaj)
                minorpairs.append(nummin)
                totalpairs.append(numpair)
                majfrac.append(nummaj / numtot)
                minfrac.append(nummin / numtot)
                totfrac.append(numpair / numtot)

            # create arrays of medians and quartiles~ 
            lower, upper = 0.5, 99.5                
            redshifts.append( redshift )
            medone.append(np.median( singleprims ))
            medtwo.append(np.median( doubleprims ))
            medtot.append(np.median( totalprims ))
            medmaj.append(np.median( majorpairs ))
            medmin.append(np.median( minorpairs ))
            medpair.append(np.median( totalpairs ))
            medmajfrac.append(np.median( majfrac ))
            medminfrac.append(np.median( minfrac ))
            medtotfrac.append(np.median( totfrac ))   
            quartone.append( np.percentile( singleprims, [lower, upper]))
            quarttwo.append( np.percentile( doubleprims, [lower, upper]))
            quarttot.append( np.percentile( totalprims, [lower, upper]))
            quartmaj.append( np.percentile( majorpairs, [lower, upper]))
            quartmin.append( np.percentile( minorpairs, [lower, upper]))
            quartpair.append( np.percentile( totalpairs, [lower, upper]))
            quartmajfrac.append( np.percentile( majfrac, [lower, upper]))
            quartminfrac.append( np.percentile( minfrac, [lower, upper]))
            quarttotfrac.append( np.percentile( totfrac, [lower, upper]))
                       
        except KeyError:
            if errorprint: print(f'skipping {snap} for KeyError. Please debug')
            continue
            
        except EmptyFile:
            if errorprint: print(f"skipping {snap}, empty file")
            continue
            
        except SkipRedshift:
            if errorprint: print(f"skipping {snap}, redshift out of range")
                
    count_dictionary = {
            "z": np.array(redshifts),
            "Median Isolated Primaries": np.array(medone),
            "Median Noniso Primaries": np.array(medtwo),
            "Median Total Primaries": np.array(medtot),
            "Median Major Pairs": np.array(medmaj),
            "Median Minor Pairs": np.array(medmin),
            "Median All Pairs": np.array(medpair),
            "Median Major Fraction": np.array(medmajfrac),
            "Median Minor Fraction": np.array(medminfrac),
            "Median Total Fraction": np.array(medtotfrac),
            "Quarts Isolated Primaries": np.array(quartone),
            "Quarts Noniso Primaries": np.array(quarttwo),
            "Quarts Total Primaries": np.array(quarttot),
            "Quarts Major Pairs": np.array(quartmaj),
            "Quarts Minor Pairs": np.array(quartmin),
            "Quarts All Pairs": np.array(quartpair),
            "Quarts Major Fraction": np.array(quartmajfrac),
            "Quarts Minor Fraction": np.array(quartminfrac),
            "Quarts Total Fraction": np.array(quarttotfrac)}
            
    return count_dictionary

def make_countdata(reals):
    #check if files exists
    filepath = f"{paths.path_plotdata}counts.hdf5"
    if not os.path.isfile(filepath):
        print("file does not exist...")
        print("creating file")
        f = h5py.File(filepath, 'w')
        print("file created successfully. adding header...")
        header_dict = {"1000 Reals - Quartile Range":"0.5-99.5%",
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
        dwarfs = get_counts("dwarf", reals)
        print("finished creating dwarf tables")
        massives = get_counts("massive", reals)
        print("finished creating massive tables")
              
        print("creating hdf5 structure")
        for size_data, size_name in zip([dwarfs,massives],["dwarf","massive"]):
            for key, val in size_data.items():
                val = np.array(val)
                dset = f.create_dataset(f'/{reals} Realizations/{size_name}/{key}', 
                                        shape=val.shape,
                                        dtype=val.dtype)
                dset[:] = val
                
        f.close()
        print("data saved~")
    

def get_smr(size, z, reals):
    zloc = np.where( zs['z'] == z)[0]
    sim = "TNG"
    snapshot = zs['ztng'][zloc][0] 

    pair_path = f"{sim}_{snapshot}_{reals}.hdf5"
    pair_data = h5py.File(f"{paths.path_pairs}{pair_path}", "r")
    
    pairs = pair_data["pairs"]["hydro"]

    pri_stell = np.array(pairs["Sub1 Stellar Mass"])
    sec_stell = np.array(pairs["Sub2 Stellar Mass"])
    seps = np.array(pairs["Separation"]) 

    # masks            
    pair_pri = get_primmask(pri_stell, size)
    pair_group = get_groupmask(np.array(pairs["Group Mass"]), size)
    pair_sepcut = seps > 10
    
    majors = (sec_stell/pri_stell > 1/4)
    minors = (sec_stell/pri_stell > 1/10) & (sec_stell/pri_stell < 1/4)

    pair_mask = pair_pri & pair_group & pair_sepcut

    major_mask = pair_mask & majors
    minor_mask = pair_mask & minors

    majors = np.array(pairs["Stellar Mass Ratio"])[major_mask]
    minors = np.array(pairs["Stellar Mass Ratio"])[minor_mask]

    key_dict = {"major":majors, "minor":minors}
    

    return key_dict     
    
def make_smrdata(reals):
    #check if files exists
    filepath = f"{paths.path_plotdata}smr.hdf5"
    if not os.path.isfile(filepath):
        print("file does not exist...")
        print("creating file")
        f = h5py.File(filepath, 'w')
        print("file created successfully. adding header...")
        header_dict = {"Simulation":"TNG100-1 (Hydro)"}


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
        
        for z in zs["z"]:
            dwarf = get_smr("dwarf", z, reals)
            massive = get_smr("massive", z, reals)
            print(f"finished z={z}")
              
            for size_data, size_name in zip([dwarf,massive],["dwarf","massive"]):
                for key, val in size_data.items():
                    val = np.array(val)
                    dset = f.create_dataset(f'/{reals} Realizations/z={z}/{size_name}/{key}', 
                                            shape=val.shape,
                                            dtype=val.dtype)
                    dset[:] = val
        print("data saved")        
        f.close()  
        
        
        
def get_pairfrac(reals, errorprint=False, redshiftcutoff=True):    
    snapshots = np.arange(0,100,1)
    snapshots = np.delete(snapshots, np.where(snapshots==48)[0])
    redcutoff = 4.2
        
    redshifts = []
    medfrac_majdw, medfrac_mindw, medfrac_majma, medfrac_minma  = [], [], [], []
    quartfrac_majdw, quartfrac_mindw, quartfrac_majma, quartfrac_minma = [], [], [], []
    medfrac_majdiff, medfrac_mindiff, quartfrac_majdiff, quartfrac_mindiff  = [], [], [], []

    for snap in snapshots:  
        frac_majdw, frac_mindw, frac_majma, frac_minma = [], [], [], []
        frac_majdiff, frac_mindiff = [], []
        
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
            
            unpair = pair_data["unpaired"]["hydro"]
            unpairStells = np.array(unpair["Sub1 Stellar Mass"])
            unpairGroups = np.array(unpair["Group Mass"])
            unpairReals = np.array(unpair['Realization'])
            
            pair = pair_data["pairs"]["hydro"]
            priStell = np.array(pair["Sub1 Stellar Mass"])
            secStell = np.array(pair["Sub2 Stellar Mass"])
            pairGroups = np.array(pair["Group Mass"])
            pairReals = np.array(pair["Realization"])
            seps = np.array(pair["Separation"]) 
            
            # subset masks for unpaired
            majors = (secStell/priStell > 1/4)
            minors = (secStell/priStell > 1/10) & (secStell/priStell < 1/4)
            allpairs = (majors + minors)
            pair_lowsep = (seps > 10)
            
            ## dwarfs
            unpair_pri_dwarf = get_primmask(unpairStells, "dwarf")
            unpair_group_dwarf = get_groupmask(unpairGroups, "dwarf")                
            
            pair_pri_dwarf = get_primmask(priStell, "dwarf")
            pair_group_dwarf = get_groupmask(pairGroups, "dwarf")
            
                # defining combined masks 
            unpair_mask_dwarf = unpair_pri_dwarf & unpair_group_dwarf
            primary_mask_dwarf = pair_pri_dwarf & pair_group_dwarf
            pair_mask_dwarf = pair_pri_dwarf & pair_group_dwarf & pair_lowsep & allpairs
            
            ## massive
            unpair_pri_massive = get_primmask(unpairStells, "massive")
            unpair_group_massive = get_groupmask(unpairGroups, "massive")                
            
            pair_pri_massive = get_primmask(priStell, "massive")
            pair_group_massive = get_groupmask(pairGroups, "massive")

                # defining combined masks 
            unpair_mask_massive = unpair_pri_massive & unpair_group_massive
            primary_mask_massive = pair_pri_massive & pair_group_massive
            pair_mask_massive = pair_pri_massive & pair_group_massive & pair_lowsep & allpairs
                                                          
            for real in np.unique(unpairReals):                  
                # make realization masks
                unpair_real = unpairReals == real
                pair_real = pairReals == real

                # make count values for single realization
                numone_dwarf = np.count_nonzero(unpair_mask_dwarf & unpair_real)
                numtwo_dwarf = np.count_nonzero(primary_mask_dwarf & pair_real)
                numtot_dwarf = numone_dwarf + numtwo_dwarf
                nummaj_dwarf = np.count_nonzero(pair_mask_dwarf & pair_real & majors)
                nummin_dwarf = np.count_nonzero(pair_mask_dwarf & pair_real & minors)
                numpair_dwarf = np.count_nonzero(pair_mask_dwarf & pair_real)
                
                # make count values for single realization
                numone_massive = np.count_nonzero(unpair_mask_massive & unpair_real)
                numtwo_massive = np.count_nonzero(primary_mask_massive & pair_real)
                numtot_massive = numone_massive + numtwo_massive
                nummaj_massive = np.count_nonzero(pair_mask_massive & pair_real & majors)
                nummin_massive = np.count_nonzero(pair_mask_massive & pair_real & minors)
                numpair_massive = np.count_nonzero(pair_mask_massive & pair_real)
                
                if (numtot_dwarf == 0) or (numtot_massive == 0):
                    continue
                    
                # collect vals for all reals
                frac_majdw.append( nummaj_dwarf/numtot_dwarf ) 
                frac_mindw.append( nummin_dwarf/numtot_dwarf ) 
                frac_majma.append( nummaj_massive/numtot_massive ) 
                frac_minma.append( nummin_massive/numtot_massive ) 
                frac_majdiff.append( (nummaj_massive/numtot_massive) - (nummaj_dwarf/numtot_dwarf) ) 
                frac_mindiff.append( (nummin_massive/numtot_massive) - (nummaj_dwarf/numtot_dwarf) ) 
                    
            # create arrays of medians and quartiles~ 
            lower, upper = 0.5, 99.5                
            redshifts.append( redshift )
            
            medfrac_majdw.append( np.median( frac_majdw ) )
            medfrac_mindw.append( np.median( frac_mindw ) )
            medfrac_majma.append( np.median( frac_majma ) )
            medfrac_minma.append( np.median( frac_minma ) )
            medfrac_majdiff.append( np.median( frac_majdiff ) )
            medfrac_mindiff.append( np.median( frac_mindiff ) )
            
            quartfrac_majdw.append( np.percentile( frac_majdw, [lower,upper] ) )
            quartfrac_mindw.append( np.percentile( frac_mindw, [lower,upper] ) )
            quartfrac_majma.append( np.percentile( frac_majma, [lower,upper] ) )
            quartfrac_minma.append( np.percentile( frac_minma, [lower,upper] ) )
            quartfrac_majdiff.append( np.percentile( frac_majdiff, [lower,upper] ) )
            quartfrac_mindiff.append( np.percentile( frac_mindiff, [lower,upper] ) )

                                   
        except KeyError:
            if errorprint: print(f'skipping {snap} for KeyError. Please debug')
            continue
            
        except EmptyFile:
            if errorprint: print(f"skipping {snap}, empty file")
            continue
            
        except SkipRedshift:
            if errorprint: print(f"skipping {snap}, redshift out of range")
                
    count_dictionary = {
            "z": np.array(redshifts),
        
            "Median Major Dwarf": np.array(medfrac_majdw),
            "Median Minor Dwarf": np.array(medfrac_mindw),
            "Median Major Massive": np.array(medfrac_majma),
            "Median Minor Massive": np.array(medfrac_minma),
            "Median Major Difference": np.array(medfrac_majdiff),
            "Median Minor Difference": np.array(medfrac_mindiff),
        
            "Quartile Major Dwarf": np.array(quartfrac_majdw),
            "Quartile Minor Dwarf": np.array(quartfrac_mindw),
            "Quartile Major Massive": np.array(quartfrac_majma),
            "Quartile Minor Massive": np.array(quartfrac_minma),
            "Quartile Major Difference": np.array(quartfrac_majdiff),
            "Quartile Minor Difference": np.array(quartfrac_mindiff) }
    
    return count_dictionary

def make_pairfracdata(reals):
    #check if files exists
    filepath = f"{paths.path_plotdata}smr.hdf5"
    if not os.path.isfile(filepath):
        print("file does not exist...")
        print("creating file")
        f = h5py.File(filepath, 'w')
        print("file created successfully. adding header...")
        header_dict = {"Simulation":"TNG100-1 (Hydro)"}


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
        
        ratios = get_pairratio()

        for key, val in ratios.items():
            val = np.array(val)
            dset = f.create_dataset(f'/{reals} Realizations/{key}', 
                                    shape=val.shape,
                                    dtype=val.dtype)
            dset[:] = val
        print("data saved")        
        f.close()  
    
    
    
    
def get_sepvel(reals, key, scaled=False, errorprint=False, redshiftcutoff=True): 
    snapshots = np.arange(0,100,1)
    snapshots = np.delete(snapshots, np.where(snapshots==48)[0])
    redcutoff = 4.2
        
    redshifts = []  
    med_majdw, med_mindw, med_majma, med_minma, med_majdiff, med_mindiff = [], [], [], [], [], []
    quart_majdw, quart_mindw, quart_majma, quart_minma, quart_majdiff, quart_mindiff = [], [], [], [], [], []
            
    for snap in snapshots:  
        try:
            pair_path = f"TNG_{snap}_{reals}.hdf5"
            pair_data = h5py.File(f"{paths.path_pairs}{pair_path}", "r")
            
            if np.size(pair_data) == 0:
                raise EmptyFile
                
            redshift = pair_data['Header'].attrs['Redshift']
            
            if redshiftcutoff & ( redshift > redcutoff) :
                raise SkipRedshift
                
            if (len(pair_data['pairs'][phys]['Group ID']) == 0):    
                raise EmptyFile
                
            pair = pair_data["pairs"][phys]
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

            lower, upper = 16,84         
            redshifts.append( redshift )
            
            med_majdw.append( np.median(real_majdw) )
            med_mindw.append( np.median(real_mindw) )
            med_majma.append( np.median(real_majma) )
            med_minma.append( np.median(real_minma) )
            med_majdiff.append( np.median(real_majdiff) )
            med_mindiff.append( np.median(real_mindiff) )
            
            quart_majdw.append( np.percentile( real_majdw, [lower,upper] ) )
            quart_mindw.append( np.percentile( real_mindw, [lower,upper] ) )
            quart_majma.append( np.percentile( real_majma, [lower,upper] ) )
            quart_minma.append( np.percentile( real_minma, [lower,upper] ) )
            quart_majdiff.append( np.percentile( real_majdiff, [lower,upper] ) )
            quart_mindiff.append( np.percentile( real_mindiff, [lower,upper] ) )

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

                        "Quartile Major Dwarf": np.array(quart_majdw),
                        "Quartile Minor Dwarf": np.array(quart_mindw),
                        "Quartile Major Massive": np.array(quart_majma),
                        "Quartile Minor Massive": np.array(quart_minma),
                        "Quartile Major Difference": np.array(quart_majdiff),
                        "Quartile Minor Difference": np.array(quart_mindiff) }
    
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
       
    


def get_sepveldist(reals, size, key, z, scaled=False):
    zloc = np.where( zs['z'] == z)[0]
    sim = "TNG"
    snapshot = zs['ztng'][zloc][0] 

    pair_path = f"{sim}_{snapshot}_{reals}.hdf5"
    pair_data = h5py.File(f"{paths.path_pairs}{pair_path}", "r")
    
    pairs = pair_data["pairs"]["hydro"]

    pri_stell = np.array(pairs["Sub1 Stellar Mass"])
    sec_stell = np.array(pairs["Sub2 Stellar Mass"])
    pairGroups = np.array(pairs["Group Mass"])
    pairGrRads = np.array(pairs["Group Radius"])
    seps = np.array(pairs["Separation"]) 
    vels = np.array(pairs["RelVel"]) 
    
    # masks            
    pair_pri = get_primmask(pri_stell, size)
    pair_group = get_groupmask(pairGroups, size)
    pair_sepcut = seps > 10
    
    pair_mask = pair_pri & pair_group & pair_sepcut
    
    majors = (sec_stell/pri_stell > 1/4)
    minors = (sec_stell/pri_stell > 1/10) & (sec_stell/pri_stell < 1/4)

    major_mask = pair_mask & majors
    minor_mask = pair_mask & minors
    
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

    majors = dat[major_mask]
    minors = dat[minor_mask]

    key_dict = {"major":majors, "minor":minors}
    

    return key_dict

def make_sepveldistdata(reals):
        #check if files exists
    filepath = f"{paths.path_plotdata}sepveldist.hdf5"
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

        print("finished creating dwarf tables")
        
        labels = ["Scaled Separation", "Separation", "Scaled Velocity", "Velocity"]         
        
        for z in zs["z"]:
            for size in ["dwarf","massive"]:
                ssep = get_sepveldist(reals, size, "Separation", z, scaled=True)
                sep = get_sepveldist(reals, size, "Separation", z, scaled=False)
                svel = get_sepveldist(reals, size, "RelVel", z, scaled=True)
                vel = get_sepveldist(reals, size, "RelVel", z, scaled=False)
              
                for datt, datt in zip([ssep,sep,svel,vel],labels):
                    for key, val in datt.items():
                        val = np.array(val)
                        dset = f.create_dataset(f'/{reals} Realizations/z={z}/{size_name}/datt/{key}', 
                                                shape=val.shape,
                                                dtype=val.dtype)
                    dset[:] = val
                    
            print(f"finished z={z}")        
                        
        f.close()
        print("data saved~")




########
########
########
# Here be ye graveyard
    


# if Plot2Make == "counts":
    
    
#     filepath = f"{paths.path_plotdata}counts.hdf5"
#     if os.path.isfile(filepath):
#         print("data file already exists")
        
#         f = h5py.File(f"{paths.path_plotdata}counts.hdf5", 'r+')
#         if f.get(f"{numreals} Realizations") is None:
#             # need to make data for real and add to data file
            
            
#         else:
#             continue
              
#     else:
#         print("data file does not already exist")
#         # make data file with realization
        
        
#     # plotting routine~
    
#     ## plot of median counts and differences for illustris and tng
#     fig, ax = plt.subplots(2, 2, figsize=(16,8), sharey=False, sharex=True, 
#                            gridspec_kw={'width_ratios': [1,1], 'height_ratios':[1,0.3],'wspace': 0.12,"hspace":0.05})

#     axd = ax[0][0]
#     axm = ax[0][1]
#     axddiff = ax[1][0]
#     axmdiff = ax[1][1]
#     axd.set(title="Dwarf", ylabel=r"N$_{\rm count}$ (thousands)")
#     axm.set(title="Massive")
#     axddiff.set(ylabel=r"N$_{\rm pair}$/N$_{\rm prim}$")



#     ##############
#     # dwarf plot #
#     ##############
#     axd.plot(dwarfs['z'],dwarfs['Median Total Primaries']/1000, color=palette["dwarf"], lw=3,label="Primaries")
#     axd.plot(dwarfs['z'],dwarfs['Median All Pairs']/1000, color=palette["dwarf"], lw=3, linestyle="dashed",label="Pairs")

#     axd.fill_between(dwarfs['z'], np.array(dwarfs['Quarts Total Primaries'])[:,0]/1000, 
#                      np.array(dwarfs['Quarts Total Primaries'])[:,1]/1000,
#                      color=palette["dwarf"],alpha=alphas["maj"])
#     axd.fill_between(dwarfs['z'], np.array(dwarfs['Quarts All Pairs'])[:,0]/1000, 
#                      np.array(dwarfs['Quarts All Pairs'])[:,1]/1000,
#                      color=palette["dwarf"],alpha=alphas["min"])

#     ################
#     # massive plot #
#     ################
#     axm.plot(massives['z'],massives['Median Total Primaries']/1000, 
#              color=palette["massive"], lw=3, label="Primaries")
#     axm.plot(massives['z'],massives['Median All Pairs']/1000, 
#              color=palette["massive"], lw=3, linestyle="dashed", label="Pairs")

#     axm.fill_between(massives['z'], np.array(massives['Quarts Total Primaries'])[:,0]/1000,
#                      np.array(massives['Quarts Total Primaries'])[:,1]/1000,
#                      color=palette["massive"],alpha=alphas["maj"])
#     axm.fill_between(massives['z'], np.array(massives['Quarts All Pairs'])[:,0]/1000, 
#                      np.array(massives['Quarts All Pairs'])[:,1]/1000,
#                      color=palette["massive"],alpha=alphas["min"])

#     ##############
#     # dwarf diff #
#     ##############
#     axddiff.plot(dwarfs['z'],dwarfs['Median Total Fraction'], color=palette["dwarf"], lw=3,label="Primaries")
#     axddiff.fill_between(dwarfs['z'], np.array(dwarfs['Quarts Total Fraction'])[:,0], np.array(dwarfs['Quarts Total Fraction'])[:,1],color=palette["dwarf"],alpha=alphas["maj"])

#     ################
#     # massive diff #
#     ################
#     axmdiff.plot(massives['z'],massives['Median Total Fraction'], color=palette["massive"], lw=3)
#     axmdiff.fill_between(massives['z'], np.array(massives['Quarts Total Fraction'])[:,0], np.array(massives['Quarts Total Fraction'])[:,1],color=palette["massive"],alpha=alphas["maj"])

#     ################
#     # Plot styling #
#     ################
#     for axx in ax[0]:
#         leg = axx.legend(loc="upper right", fontsize=16)

#     axm.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))

#     fig.supxlabel("Redshift")
#     # plt.savefig(f"../plots/paper1/counts_1000.png",bbox_inches='tight',facecolor="white")
#     plt.show()

    