""" Container class for encapsulating data about aerosol size distributions.

"""
import numpy as np

from .distributions import BaseDistribution, Lognorm, MultiModeLognorm


def dist_to_conc(dist, r_min, r_max, rule="trapezoid"):
    """Converts a swath of a size distribution function to an actual number
    concentration.

    Aerosol size distributions are typically reported by normalizing the
    number density by the size of the aerosol. However, it's sometimes more
    convenient to simply have a histogram of representing several aerosol
    size ranges (bins) and the actual number concentration one should expect
    in those bins. To accomplish this, one only needs to integrate the size
    distribution function over the range spanned by the bin.

    Parameters
    ----------
    dist : object implementing a ``pdf()`` method
        the representation of the size distribution
    r_min, r_max : float
        the lower and upper bounds of the size bin, in the native units of ``dist``
    rule : {'trapezoid', 'simpson', 'other'} (default='trapezoid')
        rule used to integrate the size distribution

    Returns
    -------
    float
        The number concentration of aerosol particles the given bin.

    Examples
    --------

    >>> dist = Lognorm(mu=0.015, sigma=1.6, N=850.0)
    >>> r_min, r_max = 0.00326456461236 0.00335634401598
    >>> dist_to_conc(dist, r_min, r_max)
    0.114256210943

    """
    pdf = dist.pdf
    width = r_max - r_min
    if rule == "trapezoid":
        return width * 0.5 * (pdf(r_max) + pdf(r_min))
    elif rule == "simpson":
        return (width / 6.0) * (
            pdf(r_max) + pdf(r_min) + 4.0 * pdf(0.5 * (r_max + r_min))
        )
    else:
        return width * pdf(0.5 * (r_max + r_min))


class AerosolSpecies(object):
    """Container class for organizing aerosol metadata.

    To allow flexibility with how aerosols are defined in the model, this class is
    meant to act as a wrapper to contain metadata about aerosols (their species
    name, etc), their chemical composition (particle mass, hygroscopicity, etc),
    and the particular size distribution chosen for the initial dry aerosol.
    Because the latter could be very diverse - for instance, it might be desired
    to have a monodisperse aerosol population, or a bin representation of a
    canonical size distribution - the core of this class is designed to take
    those representations and homogenize them for use in the model.

    To construct an :class:`AerosolSpecies`, only the metadata (``species`` and
    ``kappa``)  and the size distribution needs to be specified. The size distribution
    (``distribution``) can be an instance of :class:`Lognorm`, as
    long as an extra parameter ``bins``, which is an integer representing how many
    bins into which the distribution should be divided, is also passed to the
    constructor. In this  case, the constructor will figure out how to slice the
    size distribution to calculate all the aerosol dry radii and their number
    concentrations. If ``r_min`` and ``r_max`` are supplied, then the size range of
    the aerosols will be bracketed; else, the supplied ``distribution`` will contain
    a shape parameter or other bounds to use.

    Alternatively, a :class:`dict` can be passed as ``distribution`` where that
    slicing has already occurred. In this case, `distribution` must have 2 keys:
    ``r_drys`` and ``Nis``. Each of the values stored to those keys should fit the
    attribute descriptors above (although they don't need to be  arrays - they can
    be any iterable.)

    Parameters
    ----------
    species : string
        Name of aerosol species.
    distribution : { LogNorm, MultiLogNorm, dict }
        Representation of aerosol size distribution.
    kappa : float
        Hygroscopicity of species.
    rho : float, optional
        Density of dry aerosol material, kg m**-3.
    mw : float, optional
        Molecular weight of dry aerosol material, kg/mol.
    bins : int
        Number of bins in discretized size distribution.

    Attributes
    ----------
    nr : float
        Number of sizes tracked for this aerosol.
    r_drys : array of floats of length ``nr``
        Dry radii of each representative size tracked for this aerosol, m.
    rs : array of floats of length ``nr + 1``
        Edges of bins in discretized aerosol distribution representation, m.
    Nis : array of floats of length ``nr``
        Number concentration of aerosol of each representative size, m**-3.
    total_N : float
        Total number concentration of aerosol in this species, cm**-3.


    Examples
    --------

    Constructing sulfate aerosol with a specified lognormal distribution -

    >>> aerosol1 = AerosolSpecies('(NH4)2SO4', Lognorm(mu=0.05, sigma=2.0, N=300.),
    ...                           bins=200, kappa=0.6)

    Constructing a monodisperse sodium chloride distribution -

    >>> aerosol2 = AerosolSpecies('NaCl', {'r_drys': [0.25, ], 'Nis': [1000.0, ]},
    ...                          kappa=0.2)

    .. warning ::

        Throws a :class:`ValueError` if an unknown type of ``distribution`` is passed
        to the constructor, or if `bins` isn't present when ``distribution`` is
        an instance of :class:`Lognorm`.

    """

    def __init__(
        self,
        species,
        distribution,
        kappa,
        rho=None,
        mw=None,
        bins=None,
        r_min=None,
        r_max=None,
    ):
        self.species = species  # Species molecular formula
        self.kappa = kappa  # Kappa hygroscopicity parameter
        self.rho = rho  # aerosol density kg/m^3
        self.mw = mw
        self.bins = bins  # Number of bins for discretizing the size distribution

        # Handle the size distribution passed to the constructor
        self.distribution = distribution
        if isinstance(distribution, dict):
            self.r_drys = np.array(distribution["r_drys"]) * 1e-6

            # Compute boundaries for bins. To do this, assume the right
            # edge of the first bin is the geometric mean of the two smallest
            # dry radii. Then, always assume that r_dry is the geometric mean
            # of a bin and use that to back out all other edges in sequence
            if len(self.r_drys) > 1:
                mid1 = np.sqrt(self.r_drys[0] * self.r_drys[1])
                lr = (self.r_drys[0] ** 2.0) / mid1
                rs = [lr, mid1]
                for r_dry in self.r_drys[1:]:
                    rs.append(r_dry**2.0 / rs[-1])
                self.rs = np.array(rs) * 1e6
            else:
                # Truly mono-disperse, so no boundaries (we don't actually need
                # them in this case anyways)
                self.rs = None
            self.Nis = np.array(distribution["Nis"])
            self.N = np.sum(self.Nis)

        elif isinstance(distribution, Lognorm):
            # Check for missing keyword argument
            if bins is None:
                raise ValueError(
                    "Need to specify `bins` argument if passing a Lognorm "
                    "distribution"
                )

            if isinstance(bins, (list, np.ndarray)):
                self.rs = bins[:]
            else:
                if not r_min:
                    lr = np.log10(distribution.mu / (10.0 * distribution.sigma))
                else:
                    lr = np.log10(r_min)
                if not r_max:
                    rr = np.log10(distribution.mu * 10.0 * distribution.sigma)
                else:
                    rr = np.log10(r_max)
                self.rs = np.logspace(lr, rr, num=bins + 1)[:]

            nbins = len(self.rs)
            mids = np.array(
                [np.sqrt(a * b) for a, b in zip(self.rs[:-1], self.rs[1:])]
            )[0:nbins]
            self.Nis = np.array(
                [
                    0.5 * (b - a) * (distribution.pdf(a) + distribution.pdf(b))
                    for a, b in zip(self.rs[:-1], self.rs[1:])
                ]
            )[0:nbins]
            self.r_drys = mids * 1e-6

        elif isinstance(distribution, MultiModeLognorm):
            if bins is None:
                raise ValueError(
                    "Need to specify `bins` argument if passing a "
                    "MultiModeLognorm distribution"
                )

            small_mu = distribution.mus[0]
            small_sigma = distribution.sigmas[0]
            big_mu = distribution.mus[-1]
            big_sigma = distribution.sigmas[-1]

            if isinstance(bins, (list, np.ndarray)):
                self.rs = bins[:]
            else:
                if not r_min:
                    lr = np.log10(small_mu / (10.0 * small_sigma))
                else:
                    lr = np.log10(r_min)
                if not r_max:
                    rr = np.log10(big_mu * 10.0 * big_sigma)
                else:
                    rr = np.log10(r_max)

                self.rs = np.logspace(lr, rr, num=bins + 1)[:]
            nbins = len(self.rs)
            mids = np.array(
                [np.sqrt(a * b) for a, b in zip(self.rs[:-1], self.rs[1:])]
            )[0:nbins]
            self.Nis = np.array(
                [
                    0.5 * (b - a) * (distribution.pdf(a) + distribution.pdf(b))
                    for a, b in zip(self.rs[:-1], self.rs[1:])
                ]
            )[0:nbins]
            self.r_drys = mids * 1e-6

        else:
            raise ValueError(
                "Could not work with size distribution of type %r" % type(distribution)
            )

        # Correct to SI units
        # Nis: cm**-3 -> m**-3
        self.total_N = np.sum(self.Nis)
        self.Nis *= 1e6
        self.nr = len(self.r_drys)

    def stats(self):
        """Compute useful statistics about this aerosol's size distribution.

        Returns
        -------
        dict
            Inherits the values from the ``distribution``, and if ``rho``
            was provided, adds some statistics about the mass and
            mass-weighted properties.

        Raises
        ------
        ValueError
            If the stored ``distribution`` does not implement a ``stats()``
            function.
        """

        if isinstance(self.distribution, BaseDistribution):
            stats_dict = self.distribution.stats()

            if self.rho:
                stats_dict["total_mass"] = stats_dict["total_volume"] * self.rho
                stats_dict["mean_mass"] = stats_dict["mean_volume"] * self.rho
                stats_dict["specific_surface_area"] = (
                    stats_dict["total_surface_area"] / stats_dict["total_mass"]
                )

            return self.rho

        else:
            raise ValueError(
                "Could not work with size distribution of type %r"
                % type(self.distribution)
            )

    def __repr__(self):
        return "%s - %r" % (self.species, self.distribution)
