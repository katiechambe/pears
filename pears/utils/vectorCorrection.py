import numpy as np

h = 0.704 # this is not a function of redshift!

################################################################
# this script takes in the positions of the primary and        #
# secondary halos and returns the corrected separation vector! #
# between them in physical units, accounting for the           #
# periodic boundary conditions of the simulation box           #
################################################################

def vectorCorrection(primaryPosition, secondaryPosition, scale):
    boxSize = 75000*scale/h # box size in kpc, since box size is 75 cMpc/h! 
    differenceVector = (primaryPosition - secondaryPosition)*scale/h # positions are in comoving kpc/h
    differenceVector = np.where( differenceVector < -boxSize/2, differenceVector + boxSize, np.where( differenceVector > boxSize/2, differenceVector - boxSize , differenceVector ) )
    return differenceVector # returns in physical kpc

