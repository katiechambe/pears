""" 
Creates abundance matched realizations for dark matter halos

Usage:
------
To get stellar mass gives max mass:
    AbundanceMatching(maxmass,redshift,#samples).stellar_mass()

Details:
--------
Follows Moster, Naab, and White (2012)
https://arxiv.org/pdf/1205.5807.pdf
equations (11-14) and 
"""

__author__ = "Katie Chamberlain"
__status__ = "Beta - forever~"
__date__   = "May 2019 - edited Oct. 2021"

from numpy.random import normal

class AbundanceMatching:
    def __init__(self, maxmass, redshift, samples):
        """
        Samples from gaussian distributed abundance matching relation 

        Parameters:
        -----------
        maxmass: float
            the dark matter maximum halo mass ever achieved
            (in Msun, NOT 1e10Msun)
        redshift:
            redshift at the present snapshot
            NOTE - not the snapshot at max mass
        samples: int
            number of realizations to sample from relationship  
        """
        
        self.maxmass = maxmass
        self.z = redshift
        self.samples = samples
        
    def logM1(self):
        """eq. 11"""
        M10      = 11.59
        M10range = 0.236
        M11      = 1.195
        M11range = 0.353

        # picking gaussian distributed variables
        M10Gauss = normal(M10, M10range,self.samples)
        M11Gauss = normal(M11, M11range, self.samples)
        return M10Gauss + M11Gauss*(self.z/(1+self.z))
    
    def N(self):
        """eq. 12"""    
        N10      = 0.0351
        N10range = 0.0058
        N11      = -0.0247
        N11range = 0.0069

        # picking gaussian distributed variables
        N10Gauss = normal(N10, N10range,self.samples)
        N11Gauss = normal(N11, N11range, self.samples)
        return N10Gauss + N11Gauss*(self.z/(1+self.z))

    def beta(self):
        """eq. 13"""
        beta10      = 1.376
        beta10range = 0.153
        beta11      = -0.826
        beta11range = 0.225

        # picking gaussian distributed variables
        beta10Gauss = normal(beta10, beta10range,self.samples),
        beta11Gauss = normal(beta11, beta11range, self.samples)
        return beta10Gauss + beta11Gauss*(self.z/(1+self.z))   

    def gamma(self):
        """eq. 14"""
        gamma10      = 0.608
        gamma10range = 0.059
        gamma11      = 0.329
        gamma11range = 0.173
        
        # picking gaussian distributed variables
        gamma10Gauss = normal(gamma10, gamma10range,self.samples)
        gamma11Gauss = normal(gamma11, gamma11range, self.samples)
        return gamma10Gauss + gamma11Gauss*(self.z/(1+self.z))

    def mass_ratio(self):
        """
        stellar to halo mass ratio

        Returns:
        --------
            Stellar mass to halo mass ratio

        eq.2 in Moster
        """
        M1 = 10**self.logM1()
        A = (self.maxmass/M1)**(-self.beta())
        B = (self.maxmass/M1)**(self.gamma())
        SHMratio = 2*self.N()*(A+B)**-1
        return SHMratio

    def stellar_mass(self):
        """ 
        stellar mass

        Returns:
        --------
            stellar mass in Msun
        """
        return self.maxmass*self.mass_ratio()

 
