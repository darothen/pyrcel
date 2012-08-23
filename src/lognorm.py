#!/usr/bin/env python 
#
# Daniel Rothenberg, 2012-8-17

"""
Class and utilities for manipulating lognormal distributions.
"""
__docformat__ = "restructuredtext"

import numpy as np
from scipy.special import erf

def fit_distributions(N_total, M_total, sigma, rho=1.0):
    """
    DEPRECATED
    
    Fits number and mass distributions given the total mass and number concetration and
    shape parameter.
    
    Need to ensure that N_total and M_total are both concentrations in the same unit volume!
    To do this, assume that rho is provided in g/cm^3. there are 1e6 micograms/micrometer^3,
    so in the computation given here, rho will be divided by 1e6 to get the proper units.
    Need this in micrometers because we want to report the final answer for diameter in 
    micrometers.
    
    :param N_total:
        The total particle number concentration, in particles/cm^3
    :param M_total:
        The total particle mass concentration, in micrograms/cm^3
    :param sigma:
        An assumed shape parameter for the distributions being fit.
    :param rho:
        The density of an individual particle, in g/cm^3. Default supplied is 1.0 (liq. water)
    
    :returns:
        A 2-tuple containing the number and mass distributions.
    
    """
    count_median_diameter = ((6.0*M_total)/(np.pi*N_total*(rho/1e6)))**(1./3.)*np.exp((- 3./2.)*sigma**2.)
    mass_median_diameter = count_median_diameter*np.exp(3.0*np.log(sigma)**2.)

    number_distribution = Lognorm(count_median_diameter, sigma, N_total)
    mass_distribution = Lognorm(mass_median_diameter, sigma, M_total)
    
    return number_distribution, mass_distribution

class Lognorm(object):
    """Lognormal distribution class.
    
    An instance of 'Lognorm' contains a construction of a lognormal distribution
    and the utilities necessary for computing statistical functions associated
    with that distribution. Instance variables are fundamental parameters which
    define specific lognormal distributions.
    
    :ivar mu:
        The *mu* (r_i, Nenes et al, 2001) parameter of a lognormal distribution, 
        otherwise the number mode radius.
    
    :ivar sigma:
        The characteristic shape parameter *sigma* of a lognormal distribution,
        otherwise the geometric standard deviation.
    
    :ivar N:
        A scaling factor for a lognormal distribution. In log-space, a lognormal 
        distribution's PDF curve should contain an area equal to 1. This optional
        factor lets one scale the area underneath the lognormal PDF curve. In most
        cases, this term represents the total aerosol number concentration.
        Default - 1.0
    """
    
    def __init__(self, mu, sigma, N=1.0):
        self.mu = mu
        self.sigma = sigma
        self.N = N
        
        # Compute moments
        self.median = self.mu
        self.mean = self.mu*np.exp(0.5*self.sigma**2)
        
    def cdf(self, x):
        erf_arg = (np.log(x/self.mu))/(np.sqrt(2.0)*np.log(self.sigma))
        return (self.N/2.0)*(1.0+erf(erf_arg))
        
    def pdf(self, x):
        """
        Probability density function
        """
        scaling = self.N/(np.sqrt(2.0*np.pi)*np.log(self.sigma))
        exponent = ((np.log(x/self.mu))**2)/(2.0*(np.log(self.sigma))**2)
        #return (scaling/x)*np.exp(-exponent)
        return (scaling/x)*np.exp(-exponent)
    
    def moment(self, k):
        geo_number_mean_volume = ((4./3.)*np.pi*self.mu**3)**k
        exponent = ((9./2.)*(k**2)*(np.log(self.sigma))**2)
        return geo_number_mean_volume*self.N*np.exp(exponent)
            
    def stats(self):
        """
        Computes relevant statistics for a lognormal distribution
        """
        stats_dict = dict()
        stats_dict['mean'] = self.mu*np.exp(0.5*self.sigma**2)
        
        
        return stats_dict
        

    def __repr__(self):
        return "Lognorm| mu = %2.2e, sigma = %2.2e, Total = %2.2e |" % (self.mu, self.sigma, self.N)