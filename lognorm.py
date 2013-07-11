"""
.. module:: parcel
    :synopsis: Class and utilities for manipulating lognormal distributions.

.. moduleauthor:: Daniel Rothenberg <darothen@mit.edu>

"""

__docformat__ = "reStructuredText"

import numpy as np
from scipy.special import erf

jaenicke_distributions = {
    "Urban": MultiModeLognorm(mus=(0.0065, 0.007, 0.025), 
                              sigmas=(10**0.245, 10**0.666, 10**0.337),
                              Ns=(9.93e4, 1.11e3, 3.64e4)),
}

class MultiModeLognorm(object):
    """Multimode lognormal distribution class.

    Container for multiple Lognorm classes representing a full aerosol size
    distribution.
    """

    def __init__(self, mus, sigmas, Ns, base="e"):
        self.mus = np.array(sorted(mus))
        self.sigmas = np.array(sorted(sigmas))
        self.Ns = np.array(sorted(Ns))
        self.base = base

        self.lognorms = []
        for mu, sigma, N in zip(self.mus, self.sigmas, self.Ns):
            mode_dist = Lognorm(mu, sigma, N, base)
            self.lognorms.append(mode_dist)

    def cdf(self, x):
        return np.sum([d.cdf(x) for d in self.lognorms])

    def pdf(self, x):
        return np.sum([d.pdf(x) for d in self.lognorms])

    def __repr__(self):
        mus_str = "("+", ".join("%2.2e" % mu for mu in self.mus)+")"
        sigmas_str = "("+", ".join("%2.2e" % sigma for sigma in self.sigmas)+")"
        Ns_str = "("+", ".join("%2.2e" % N for N in self.Ns)+")"
        return "MultiModeLognorm| mus = %s, sigmas = %s, Totals = %s |" % (mus_str, sigmas_str, Ns_str)

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

    def __init__(self, mu, sigma, N=1.0, base="e"):
        self.mu = mu
        self.sigma = sigma
        self.N = N

        self.base = base
        if self.base in ["e", np.e]:
            self.log = np.log
        elif self.base == 10.:
            self.log = np.log10
        else:
            self.log_base = np.log(self.base)
            self.log = lambda x: np.log(x)/self.log_base

        # Compute moments
        self.median = self.mu
        self.mean = self.mu*np.exp(0.5*self.sigma**2)

    def cdf(self, x):
        erf_arg = (self.log(x/self.mu))/(np.sqrt(2.0)*self.log(self.sigma))
        return (self.N/2.0)*(1.0+erf(erf_arg))

    def pdf(self, x):
        """
        Probability density function
        """
        scaling = self.N/(np.sqrt(2.0*np.pi)*self.log(self.sigma))
        exponent = ((self.log(x/self.mu))**2)/(2.0*(self.log(self.sigma))**2)
        #return (scaling/x)*np.exp(-exponent)
        return (scaling/x)*np.exp(-exponent)

    def moment(self, k):
        geo_number_mean_volume = ((4./3.)*np.pi*self.mu**3)**k
        exponent = ((9./2.)*(k**2)*(self.log(self.sigma))**2)
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
