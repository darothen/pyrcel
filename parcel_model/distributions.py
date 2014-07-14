"""
.. module:: parcel
    :synopsis: Class and utilities for manipulating lognormal distributions.

.. moduleauthor:: Daniel Rothenberg <darothen@mit.edu>

"""

__docformat__ = "reStructuredText"

import numpy as np
from scipy.special import erf

class MultiModeLognorm(object):
    """Multimode lognormal distribution class.

    Container for multiple Lognorm classes representing a full aerosol size
    distribution.
    """

    def __init__(self, mus, sigmas, Ns, base="e"):

        dist_params = zip(mus, sigmas, Ns)
        from operator import itemgetter
        dist_params = sorted(dist_params, key=itemgetter(0))

        self.mus, self.sigmas, self.Ns = zip(*dist_params)

        self.base = base

        self.lognorms = []
        for mu, sigma, N in zip(self.mus, self.sigmas, self.Ns):
            mode_dist = Lognorm(mu, sigma, N, base)
            self.lognorms.append(mode_dist)

    def cdf(self, x):
        return np.sum([d.cdf(x) for d in self.lognorms], axis=0)

    def pdf(self, x):
        return np.sum([d.pdf(x) for d in self.lognorms], axis=0)

    def __repr__(self):
        mus_str = "("+", ".join("%2.2e" % mu for mu in self.mus)+")"
        sigmas_str = "("+", ".join("%2.2e" % sigma for sigma in self.sigmas)+")"
        Ns_str = "("+", ".join("%2.2e" % N for N in self.Ns)+")"
        return "MultiModeLognorm| mus = %s, sigmas = %s, Totals = %s |" % \
               (mus_str, sigmas_str, Ns_str)


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

## Source = Aerosol-Cloud-Climate Interactions by Peter V. Hobbs, pg. 14
jaenicke_distributions = {
    "Polar": MultiModeLognorm(mus=(0.0689, 0.375, 4.29),
                              sigmas=(10**0.245, 10**0.300, 10**0.291),
                              Ns=(21.7, 0.186, 3.04e-4), base=10.),
    "Urban": MultiModeLognorm(mus=(0.00651, 0.00714, 0.0248),
                              sigmas=(10.**0.245, 10.**0.666, 10.**0.337),
                              Ns=(9.93e4, 1.11e3, 3.64e4), base=10.),
    "Background": MultiModeLognorm(mus=(0.0036, 0.127, 0.259),
                                   sigmas=(10.**0.645, 10.**0.253, 10.**0.425),
                                   Ns=(129., 59.7, 63.5), base=10.),
    "Maritime": MultiModeLognorm(mus=(0.0039, 0.133, 0.29),
                                 sigmas=(10.**0.657, 10.**0.210, 10.**0.396),
                                 Ns=(133., 66.6, 3.06), base=10.),
    "Remote Continental": MultiModeLognorm(mus=(0.01, 0.058, 0.9),
                                           sigmas=(10.**0.161, 10.**0.217, 10.**0.38),
                                           Ns=(3.2e3, 2.9e3, 0.3), base=10.),
    "Rural": MultiModeLognorm(mus=(0.00739, 0.0269, 0.0149),
                              sigmas=(10.**0.225, 10.**0.557, 10.**0.266),
                              Ns=(6.65e3, 147., 1990.), base=10.),
}