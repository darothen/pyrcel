""" Collection of classes for representing aerosol size distributions.

Most commonly, one would use the :class:`Lognorm` distribution. However,
for the sake of completeness, other canonical distributions will be
included here, with the notion that this package could be extended to
describe droplet size distributions or other collections of objects.

"""
from abc import ABCMeta, abstractmethod

import numpy as np
from scipy.special import erf, erfinv


class BaseDistribution(metaclass=ABCMeta):
    """Interface for distributions, to ensure that they contain a pdf method."""

    @abstractmethod
    def cdf(self, x):
        """Cumulative density function"""

    @abstractmethod
    def pdf(self, x):
        """Distribution density function."""

    @property
    @abstractmethod
    def stats(self):
        pass

    @abstractmethod
    def __repr__(self):
        """Representation function."""


class Gamma(BaseDistribution):
    """Gamma size distribution"""

    # TODO: Implement Gamma size distribution
    pass


class Lognorm(BaseDistribution):
    """Lognormal size distribution.

    An instance of :class:`Lognorm` contains a construction of a lognormal distribution
    and the utilities necessary for computing statistical functions associated
    with that distribution. The parameters of the constructor are invariant with respect
    to what length and concentration unit you choose; that is, if you use meters for
    ``mu`` and cm**-3 for ``N``, then you should keep these in mind when evaluating
    the :func:`pdf` and :func:`cdf` functions and when interpreting moments.

    Parameters
    ----------
    mu : float
        Median/geometric mean radius, length unit.
    sigma : float
        Geometric standard deviation, unitless.
    N : float, optional (default=1.0)
        Total number concentration, concentration unit.

    Attributes
    ----------
    median, mean : float
        Pre-computed statistical quantities

    Methods
    -------
    pdf(x)
        Evaluate size distribution density dN/dx (not logarithmic) at a particular radius x
    pdfloge(x)
        Evaluate size distribution logarithmic density dN/dln(x) at a particular radius x
    pdflog10(x)
        Evaluate size distribution logarithmic density dN/dlog(x) at a particular radius x
    cdf(x)
        Evaluate cumulative concentration up to a particular radius x
    moment(k)
        Compute the *k*-th moment of the lognormal distribution.

    """

    def __init__(self, mu, sigma, N=1.0, base=np.e):
        self.mu = mu
        self.sigma = sigma
        self.N = N

        self.base = np.e
        self.log = np.log

        # Compute moments
        self.median = self.mu
        self.mean = self.mu * np.exp(0.5 * self.sigma**2)

    def invcdf(self, y):
        """Inverse of cumulative density function.

        Parameters
        ----------
        y : float
            CDF value, between (0, 1)

        Returns
        -------
        radius (same unit as mu) corresponding to given CDF value

        """

        if (np.any(y) < 0) or (np.any(y) > 1):
            raise ValueError("y must be between (0, 1)")

        erfinv_arg = 2.0 * y / self.N - 1.0
        return self.mu * np.exp(np.log(self.sigma) * np.sqrt(2.0) * erfinv(erfinv_arg))

    def cdf(self, x):
        """Cumulative density function

        .. math::
            \\text{CDF} = \\frac{N}{2}\\left(1.0 + \\text{erf}(\\frac{\ln(x/\mu)}{\sqrt{2}\ln{\sigma}}) \\right)

        Parameters
        ----------
        x : float
            Radius (must match unit of mu)

        Returns
        -------
        cumulative concentration up to radius x

        """
        erf_arg = (np.log(x / self.mu)) / (np.sqrt(2.0) * np.log(self.sigma))
        return (self.N / 2.0) * (1.0 + erf(erf_arg))

    def pdf(self, x):
        """Distribution density function dN/dx (not logarithmic)

        .. math::
            \\text{pdf} = \\frac{N}{\sqrt{2\pi}\ln(\sigma) x}\exp\\left( -\\frac{\ln^2(x/\mu)}{2\ln^2\sigma} \\right)

        Parameters
        ----------
        x : float
            Radius (must match unit of mu)

        Returns
        -------
        distribution density value at radius x

        """
        scaling = self.N / (np.sqrt(2.0 * np.pi) * np.log(self.sigma))
        exponent = ((np.log(x / self.mu)) ** 2) / (2.0 * (np.log(self.sigma)) ** 2)
        return (scaling / x) * np.exp(-exponent)

    def pdfloge(self, x):
        """Distribution density function dN/dln(x) (natural logarithm)

        .. math::
            \\text{pdf} = \\frac{N}{\sqrt{2\pi}\ln\sigma}\exp\\left( -\\frac{\ln^2(x/\mu)}{2\ln^2\sigma} \\right)

        Parameters
        ----------
        x : float
            Radius (must match unit of mu)

        Returns
        -------
        distribution logarithmic density value at radius x

        """
        scaling = self.N / (np.sqrt(2.0 * np.pi) * np.log(self.sigma))
        exponent = ((np.log(x / self.mu)) ** 2) / (2.0 * (np.log(self.sigma)) ** 2)
        return scaling * np.exp(-exponent)

    def pdflog10(self, x):
        """Distribution density function dN/dlog(x) (base 10 logarithm)

        .. math::
            \\text{pdf} = \\frac{N\ln 10}{\sqrt{2\pi}\ln\sigma}\exp\\left( -\\frac{\ln^2(x/\mu)}{2\ln^2\sigma} \\right)

        Parameters
        ----------
        x : float
            Radius (must match unit of mu)

        Returns
        -------
        distribution logarithmic density value at radius x

        """
        scaling = self.N / (np.sqrt(2.0 * np.pi) * np.log(self.sigma))
        exponent = ((np.log(x / self.mu)) ** 2) / (2.0 * (np.log(self.sigma)) ** 2)
        return np.log(10) * scaling * np.exp(-exponent)

    def moment(self, k):
        """Compute the k-th moment of the lognormal distribution

        .. math::
            F(k) = N\mu^k\exp\\left( \\frac{k^2}{2} \ln^2 \sigma \\right)

        Parameters
        ----------
        k : int
            Moment to evaluate

        Returns
        -------
        moment of distribution

        """
        scaling = (self.mu**k) * self.N
        exponent = ((k**2) / 2.0) * (np.log(self.sigma)) ** 2
        return scaling * np.exp(exponent)

    def stats(self):
        """Compute useful statistics for a lognormal distribution

        Returns
        -------
        dict
            Dictionary containing the stats ``mean_radius``, ``total_diameter``,
            ``total_surface_area``, ``total_volume``, ``mean_surface_area``,
            ``mean_volume``, and ``effective_radius``

        """
        stats_dict = dict()
        stats_dict["mean_radius"] = self.mu * np.exp(0.5 * (np.log(self.sigma))**2)

        stats_dict["total_diameter"] = self.N * stats_dict["mean_radius"]
        stats_dict["total_surface_area"] = 4.0 * np.pi * self.moment(2.0)
        stats_dict["total_volume"] = (4.0 * np.pi / 3.0) * self.moment(3.0)

        stats_dict["mean_surface_area"] = stats_dict["total_surface_area"] / self.N
        stats_dict["mean_volume"] = stats_dict["total_volume"] / self.N

        stats_dict["effective_radius"] = (
            stats_dict["total_volume"] / stats_dict["total_surface_area"]
        )

        return stats_dict

    def __repr__(self):
        return "Lognorm | mu = {:2.2e}, sigma = {:2.2e}, Total = {:2.2e} |".format(
            self.mu, self.sigma, self.N
        )


class MultiModeLognorm(BaseDistribution):
    """Multimode lognormal distribution class.

    Container for multiple Lognorm classes representing a full aerosol size
    distribution.
    """

    def __init__(self, mus, sigmas, Ns, base=np.e):
        dist_params = list(zip(mus, sigmas, Ns))
        from operator import itemgetter

        dist_params = sorted(dist_params, key=itemgetter(0))

        self.mus, self.sigmas, self.Ns = list(zip(*dist_params))

        self.base = np.e

        self.lognorms = []
        for mu, sigma, N in zip(self.mus, self.sigmas, self.Ns):
            mode_dist = Lognorm(mu, sigma, N)
            self.lognorms.append(mode_dist)

    def cdf(self, x):
        return np.sum([d.cdf(x) for d in self.lognorms], axis=0)

    def pdf(self, x):
        return np.sum([d.pdf(x) for d in self.lognorms], axis=0)

    def stats(self):
        """Compute useful statistics for a multi-mode lognormal distribution

        TODO: Implement multi-mode lognorm stats
        """
        raise NotImplementedError()

    def __repr__(self):
        mus_str = "(" + ", ".join("%2.2e" % mu for mu in self.mus) + ")"
        sigmas_str = "(" + ", ".join("%2.2e" % sigma for sigma in self.sigmas) + ")"
        Ns_str = "(" + ", ".join("%2.2e" % N for N in self.Ns) + ")"
        return "MultiModeLognorm| mus = {}, sigmas = {}, Totals = {} |".format(
            mus_str, sigmas_str, Ns_str
        )


#: Single mode aerosols
# TODO: Re-factor saved size distributions into a resource like 'constants'

FN2005_single_modes = {
    "SM1": Lognorm(0.025 / 2, 1.3, 100.0),
    "SM2": Lognorm(0.025 / 2, 1.3, 500.0),
    "SM3": Lognorm(0.05 / 2, 1.8, 500.0),
    "SM4": Lognorm(0.25 / 2, 1.8, 100.0),
    "SM5": Lognorm(0.75 / 2, 1.8, 1000.0),
}

NS2003_single_modes = {
    "SM1": Lognorm(0.02 / 2, 2.5, 200.0),
    "SM2": Lognorm(0.02 / 2, 2.5, 1000.0),
    "SM3": Lognorm(0.02 / 2, 1.5, 1000.0),
    "SM4": Lognorm(0.2 / 2, 2.5, 200.0),
    "SM5": Lognorm(0.02 / 2, 2.5, 10000.0),
}

whitby_distributions = {
    # name: [nucleation, accumulation, coarse]
    #        mu = micron, N = cm**-3
    "marine": [
        Lognorm(0.01 / 2.0, 1.6, 340.0),
        Lognorm(0.07 / 2.0, 2.0, 60.0),
        Lognorm(0.62 / 2.0, 2.7, 3.1),
    ],
    "continental": [
        Lognorm(0.016 / 2.0, 1.6, 1000.0),
        Lognorm(0.068 / 2.0, 2.1, 800.0),
        Lognorm(0.92 / 2.0, 2.2, 0.72),
    ],
    "background": [
        Lognorm(0.01 / 2.0, 1.7, 6400.0),
        Lognorm(0.076 / 2.0, 2.0, 2300.0),
        Lognorm(1.02 / 2.0, 2.16, 3.2),
    ],
    "urban": [
        Lognorm(0.014 / 2.0, 1.8, 10600.0),
        Lognorm(0.054 / 2.0, 2.16, 32000.0),
        Lognorm(0.86 / 2.0, 2.21, 5.4),
    ],
}

# Source = Aerosol-Cloud-Climate Interactions by Peter V. Hobbs, pg. 14
jaenicke_distributions = {
    "Polar": MultiModeLognorm(
        mus=(0.0689, 0.375, 4.29),
        sigmas=(10**0.245, 10**0.300, 10**0.291),
        Ns=(21.7, 0.186, 3.04e-4),
    ),
    "Urban": MultiModeLognorm(
        mus=(0.00651, 0.00714, 0.0248),
        sigmas=(10.0**0.245, 10.0**0.666, 10.0**0.337),
        Ns=(9.93e4, 1.11e3, 3.64e4),
    ),
    "Background": MultiModeLognorm(
        mus=(0.0036, 0.127, 0.259),
        sigmas=(10.0**0.645, 10.0**0.253, 10.0**0.425),
        Ns=(129.0, 59.7, 63.5),
    ),
    "Maritime": MultiModeLognorm(
        mus=(0.0039, 0.133, 0.29),
        sigmas=(10.0**0.657, 10.0**0.210, 10.0**0.396),
        Ns=(133.0, 66.6, 3.06),
    ),
    "Remote Continental": MultiModeLognorm(
        mus=(0.01, 0.058, 0.9),
        sigmas=(10.0**0.161, 10.0**0.217, 10.0**0.38),
        Ns=(3.2e3, 2.9e3, 0.3),
    ),
    "Rural": MultiModeLognorm(
        mus=(0.00739, 0.0269, 0.0419),
        sigmas=(10.0**0.225, 10.0**0.557, 10.0**0.266),
        Ns=(6.65e3, 147.0, 1990.0),
    ),
}
