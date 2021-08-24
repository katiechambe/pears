import numpy as np

################################################################
# this script takes in the position of a subhalo and returns   #
# the coordinates of a box centered around the subhalo with    #
# 10 ckpc/h sides and accounts for the periodic boundary       #
# conditions of the simulation box                             #
################################################################

def subboxCoordinates(position):
    '''
    inputs: array of x,y,z, position of subhalos in ckpc/h, scale of snapshot
    --- 
    outputs: xLo, xHi, yLo, yHi, zLo, zHi in ckpc/h
    '''
    boxSize = 75000 # since box size is 75000 comoving kpc/h
    subboxSize = 10000 # in comoving kpc/h
    
    #xPos, yPos, zPos = position[:,0], position[:,1], position[:,2]
    xPos, yPos, zPos = position[0], position[1], position[2]

    # defining subbox coordinates 
    xLo = xPos - subboxSize
    xHi = xPos + subboxSize

    yLo = yPos - subboxSize
    yHi = yPos + subboxSize

    zLo = zPos - subboxSize
    zHi = zPos + subboxSize

    #correct for periodic box edges
    xLo = np.where(xLo < 0, xLo + boxSize, xLo)
    xHi = np.where(xHi > boxSize, xHi - boxSize, xHi)

    yLo = np.where(yLo < 0, yLo + boxSize, yLo)
    yHi = np.where(yHi > boxSize, yHi - boxSize, yHi)

    zLo = np.where(zLo < 0, zLo + boxSize, zLo)
    zHi = np.where(zHi > boxSize, zHi - boxSize, zHi)

    return xLo, xHi, yLo, yHi, zLo, zHi # returns in ckpc/h!! 
