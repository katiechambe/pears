import numpy as np

'''
vectorCorrection corrects the separation vector between two subhalos
accounting for the periodic boundary conditions of the simulation box. 
Returns the corrected separation vector between the subhalos
Note: 
----- 
Input all values in physical units (kpc)
'''

def vectorCorrection(primaryPosition, secondaryPosition, boxsize):
    '''
    Note:
    -----
    All values should be input in physical units! including box size
    '''
    differenceVector = (primaryPosition - secondaryPosition) 

    corrected = np.where(differenceVector < -boxsize/2, 
                                differenceVector + boxsize, 
                                np.where(differenceVector > boxsize/2, 
                                         differenceVector - boxsize, 
                                         differenceVector
                                         ) 
                                )
    return corrected 

