#!/usr/bin/env python 
#
# Daniel Rothenberg, 2012-8-17

"""
Class and utilities for manipulating lognormal distributions.
"""
__docformat__ = "restructuredtext"

import numpy as np
from scipy.special import erf

class MultiModeLognorm:
    """Multimode lognormal distribution class.
    
    Container for multiple Lognorm classes representing a full aerosol size
    distribution.    
    """
    
    def __init__(self, mus, sigmas, Ns):
        self.mus = np.array(mus)
        self.sigmas = np.array(sigmas)
        self.Ns = np.array(Ns)
        
        self.lognorms = []
        for mu, sigma, N in zip(self.mus, self.sigmas, self.Ns):
            mode_dist = Lognorm(mu, sigma, N)
            self.lognorms.append(mode_dist)
            
    def cdf(self, x):
        return np.sum([d.cdf(x) for d in self.lognorms])

    def pdf(self, x):
        return np.sum([d.pdf(x) for d in self.lognorms])

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