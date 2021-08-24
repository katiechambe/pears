
##################################################################
# This script creates a function that calculates the hill radius #
# for a galaxy of mass m at some distance r away from another    #
# galaxy of mass M.                                              #
# -------------------------------------------------------------- #
# Katie Chamberlain - Dec 2020                                   #
##################################################################

def hillRadius(mass, separation, massesExt):
    '''
    from: arXiv:0803.4211 Hahn+ 2009 eq. 10
    ---
    inputs: 
    - mass of group under consideration (mass)
    - separation to external groups (separation, kpc)
    - masses of external groups (massesExt)
    returns: 
    - the hill radius in physical kpc
    ''' 
    hillRadius = separation * ( mass / (3* massesExt) )**(1/3)

    return hillRadius
