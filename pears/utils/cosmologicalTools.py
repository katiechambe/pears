"""
Cosmological tools
"""

import numpy as np
from scipy.integrate import quad
import matplotlib.pyplot as plt


class CosmologicalTools:
    def __init__(self, OmegaM=0.2726, OmegaL=0, OmegaR=0, h=0.704):

        # initializing the cosmology
        self.OmegaM = OmegaM  # matter density parameter
        self.OmegaL = OmegaL  # dark energy density parameter
        self.OmegaR = OmegaR  # radiation density parameter
        self.OmegaK = 1.0 - (OmegaM + OmegaL + OmegaR)  # curvature parameter
        self.h = h  # normalization for the hubble parameter
        self.H0 = h * 100  # hubble constant at z=0 100*h km/s/Mpc

        # physical constants
        self.c = 299792.458  # km/s
        self.G = 6.67408e-11  # m^3/(kg s^2)

        # to help determine distance measure
        if self.OmegaK > 0:
            k = -1
            self.Rc = np.sqrt((-k * self.c ** 2) / ((self.H0 ** 2) * self.OmegaK))  # radius of curvature

    ###									    ###
    ###		Here are cosmography tools		###
    ###										###

    def HubbleParameterZ(self, z):
        """
        Hubble parameter as a function of redshift
        Redshift can be entered as a number or an array
        Returns in units of km/s/Mpc
        """
        Omz = self.OmegaM * (1 + z) ** 3
        Olz = self.OmegaL
        Orz = self.OmegaR * (1 + z) ** 4
        Okz = self.OmegaK * (1 + z) ** 2
        Hz = self.H0 * np.sqrt(Omz + Olz + Orz + Okz)
        return Hz

    def OmegaMZ(self, z):
        """
        Matter density parameter as a function of z
        Redshift can be entered as a number or an array
        """
        H = self.H0 / self.HubbleParameterZ(z)
        Omz = self.OmegaM * (1 + z) ** 3
        return Omz * (H ** 2)

    def OmegaLZ(self, z):
        """
        Dark energy density parameter as a function of z
        Redshift can be entered as a number or an array
        """
        H = self.H0 / self.HubbleParameterZ(z)
        Olz = self.OmegaL
        return Olz * (H ** 2)

    def OmegaRZ(self, z):
        """
        Dark energy density parameter as a function of z
        Redshift can be entered as a number or an array
        """
        H = self.H0 / self.HubbleParameterZ(z)
        Orz = self.OmegaR * (1 + z) ** 4
        return Orz * (H ** 2)

    def OmegaKZ(self, z):
        """
        Curvature parameter as a function of z
        Redshift can be entered as a number or an array
        """
        H = self.H0 / self.HubbleParameterZ(z)
        Okz = self.OmegaK * (1 + z) ** 2
        return Okz * (H ** 2)

    def comovingDistance(self, z):
        """
        Calculates the comoving distance as a fucntion of redshift
        Redshift must be a number (not an array!)
        Returns the comoving distance between 0 and z in Mpc
        """

        def integrand(x):
            return self.c / self.HubbleParameterZ(x)

        return quad(integrand, 0, z)[0]

    def distanceMeasure(self, z):
        """
        Calculates the distance measure in the case that the universe is open or flat
        Redshift must be a number
        Returns in Mpc
        """
        if self.OmegaK > 0:
            return self.Rc * np.sinh(self.comovingDistance(z) / self.Rc)
        else:
            return self.comovingDistance(z)

    def angularDiameter(self, z):
        """
        Angular diameter distance as a function of redshift
        Redshift must be a number - integrates between 0 and z
        Returns in Mpc/rad
        """
        return self.distanceMeasure(z) / (1 + z)

    def luminosityDistance(self, z):
        """
        Luminosity distance as a function of redshift
        Redshift must be a number
        Returns in Mpc
        """
        return self.distanceMeasure(z) * (1 + z)

    def distanceModulus(self, z):
        """
        Distance modulus as z
        Redshift must be a number
        """
        return 5 * np.log10(self.luminosityDistance(z) * 1e6 / 10)  # convert to pc in argument

    def lookbackTime(self, z):
        # Lookback time at z
        # Redshift must be a number
        # From z=0 to z and returns in gigayears

        def integrand(x):
            return (self.HubbleParameterZ(x) * (1 + x)) ** (-1)

        return quad(integrand, 0, z)[0] * 9.77799e2

    ###								###
    ###		Here are halo tools		###
    ###								###

    def deltavir(self, z):
        # ---
        # Input(s): redshift (might not be necessary?)
        # Return(s): virial overdensity (dimensionless)
        # ...
        # From: Bryan and Norman 1998 arXiv:astro-ph/9710107v1
        # ---
        x = self.OmegaMZ(z) - 1
        A = 18 * np.pi ** 2
        B = 82 * x
        C = 39 * x ** 2
        deltac = A + B - C
        return deltac / self.OmegaMZ(z)

    def cvir(self, mvir):
        """
        Input(s): mvir in Msun
        Return(s): virial concentration (dimensionless)
        ...
        From: Klypin 2011 arxiv:1002.3660v4
        """
        # TODO: this should be redshift dependent...

        return 9.6 * ((mvir * self.h / 1e12)**(-0.075))
        
    def f(self,x):
        return np.log(1 + x) - (x / (1 + x))

    def c200(self, cvir, z):
        """
        ...
        From: van der Marel 2012 arXiv:1205.6864v1
        """
        q = 2.058  # 200/(self.OmegaMZ(z)*self.deltavir(z)) # does omega M need to be z-dependent?
        error = 1
        i = 0
        guess = cvir
        c200val = 0
        while (i < 100) or (error < 1e-10):  # take 100 steps with tolerance of 1e-10
            c200val = cvir * (self.f(guess) / (q * self.f(cvir))) ** (1. / 3.)
            error = np.abs(cvir - c200val)
            guess = c200val
            i = i + 1
        return c200val

    def m200frac(self, cvir, z):
        """
        ...
        From: van der Marel 2012 arXiv:1205.6864v1
        """
        A = self.f(self.c200(cvir, z))
        B = self.f(cvir)
        return A / B
        
    def afrac(self,c): #a/rs
        A = (2*self.f(c))**(-1/2)
        B = 1/c
        return (A-B)**(-1)
        
        
    def MHfrac(self,c): # Mh/Mvir
        A = self.afrac(c)**2
        B = 2*self.f(10)
        return A/B
        
        
        
    def rvir(self, z, mvir):
        '''
        In physical units! kpc
        ...
        From van der Marel 2012 arXiv:1205.6864v1
        '''
        # be sure mvir is in Msun, not 10^10 Msun
        A = self.deltavir(z)*self.OmegaMZ(z)/97.2
        B = mvir*self.h/1e12 # divide by 1e12 since Mvir is in Msun
        # divide by 100 if Mvir is in 10^10 Msun
        return (206/self.h)*(B/A)**(1/3)
