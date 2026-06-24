"""Container class for encapsulating data about aerosol size distributions."""

from __future__ import annotations

from typing import Any

import numpy as np
from numpy.typing import NDArray

from .distributions import BaseDistribution, Lognorm, MultiModeLognorm


def dist_to_conc(dist: BaseDistribution, r_min: float, r_max: float, rule: str = "trapezoid") -> float:
    """Convert a size distribution over a bin interval to a number concentration.

    Aerosol size distributions are typically reported as dN/dr (number density
    per unit radius). This function integrates the pdf over ``[r_min, r_max]``
    using a simple quadrature rule to obtain the actual number concentration
    in that bin.

    Parameters
    ----------
    dist : object with a ``pdf(r)`` method
        Representation of the size distribution.
    r_min, r_max : float
        Lower and upper bounds of the size bin, in the native units of ``dist``.
    rule : {'trapezoid', 'simpson', 'midpoint'}, default 'trapezoid'
        Quadrature rule.

    Returns
    -------
    float
        Number concentration of aerosol particles in the given bin.

    Examples
    --------
    >>> dist = Lognorm(mu=0.015, sigma=1.6, N=850.0)
    >>> r_min, r_max = 0.00326456461236, 0.00335634401598
    >>> dist_to_conc(dist, r_min, r_max)
    0.114256210943
    """
    pdf = dist.pdf
    width = r_max - r_min
    if rule == "trapezoid":
        return width * 0.5 * (pdf(r_max) + pdf(r_min))
    elif rule == "simpson":
        return (width / 6.0) * (pdf(r_max) + pdf(r_min) + 4.0 * pdf(0.5 * (r_max + r_min)))
    else:
        return width * pdf(0.5 * (r_max + r_min))


def _discretize_dist(dist, bins, r_min, r_max, mu_lo, sigma_lo, mu_hi, sigma_hi):
    """Discretize a continuous size distribution into bin edges, mid-radii, and Ni.

    Parameters
    ----------
    dist : BaseDistribution
        The distribution to discretize (must implement ``pdf``).
    bins : int or array-like
        If an array, used directly as bin edges (microns). If an int, that many
        logarithmically-spaced edges are generated between ``r_min``/``r_max``
        or the distribution-derived defaults.
    r_min, r_max : float or None
        Override for the lower/upper radius bounds (microns). ``None`` means
        derive from the distribution parameters.
    mu_lo, sigma_lo : float
        Median and geometric std-dev of the smallest mode, used when ``r_min``
        is not supplied.
    mu_hi, sigma_hi : float
        Median and geometric std-dev of the largest mode, used when ``r_max``
        is not supplied.

    Returns
    -------
    rs : ndarray
        Bin edges (microns), length ``n_bins + 1``.
    r_drys : ndarray
        Geometric mid-radii of each bin (metres), length ``n_bins``.
    Nis : ndarray
        Number concentration per bin (cm⁻³), length ``n_bins``.
    """
    if isinstance(bins, (list, np.ndarray)):
        rs = np.asarray(bins, dtype=float)
    else:
        log_r_min = np.log10(r_min) if r_min else np.log10(mu_lo / (10.0 * sigma_lo))
        log_r_max = np.log10(r_max) if r_max else np.log10(mu_hi * 10.0 * sigma_hi)
        rs = np.logspace(log_r_min, log_r_max, num=bins + 1)

    mids = np.sqrt(rs[:-1] * rs[1:])
    Nis = np.array([0.5 * (b - a) * (dist.pdf(a) + dist.pdf(b)) for a, b in zip(rs[:-1], rs[1:])])
    r_drys = mids * 1e-6
    return rs, r_drys, Nis


class AerosolSpecies:
    """Container for aerosol metadata and discretized size distribution.

    Wraps a species name, hygroscopicity, and size distribution into a single
    object suitable for use as a parcel-model input. The constructor accepts
    three kinds of size distributions:

    * :class:`~pyrcel.distributions.Lognorm` — a single log-normal mode; ``bins``
      must also be supplied.
    * :class:`~pyrcel.distributions.MultiModeLognorm` — multiple log-normal modes;
      ``bins`` must also be supplied.
    * ``dict`` with keys ``"r_drys"`` (microns) and ``"Nis"`` (cm⁻³) — a
      pre-discretized distribution.

    Parameters
    ----------
    species : str
        Name or molecular formula of the aerosol species.
    distribution : Lognorm, MultiModeLognorm, or dict
        Size distribution. If a ``dict``, must have keys ``"r_drys"`` and
        ``"Nis"``. If a continuous distribution, ``bins`` is required.
    kappa : float
        Kappa hygroscopicity parameter (dimensionless).
    rho : float, optional
        Density of the dry aerosol material (kg m⁻³).
    mw : float, optional
        Molecular weight of the dry aerosol material (kg mol⁻¹).
    bins : int or array-like, optional
        Number of bins (int), or explicit bin-edge array in microns.
        Required when ``distribution`` is a :class:`Lognorm` or
        :class:`MultiModeLognorm`. If an array, those edges are used directly.
    r_min, r_max : float, optional
        Override the auto-derived lower/upper radius bounds (microns) when
        computing log-spaced bins.

    Attributes
    ----------
    nr : int
        Number of size bins.
    r_drys : ndarray, shape (nr,)
        Dry radius of each representative bin (m).
    rs : ndarray or None
        Bin edges (microns) for the discretized distribution; ``None`` for a
        single-bin monodisperse input.
    Nis : ndarray, shape (nr,)
        Number concentration of each bin (m⁻³).
    total_N : float
        Total number concentration (cm⁻³).

    Raises
    ------
    ValueError
        If ``distribution`` is an unrecognised type, or if ``bins`` is missing
        for a continuous distribution.

    Examples
    --------
    Log-normal sulfate aerosol:

    >>> aerosol1 = AerosolSpecies('(NH4)2SO4', Lognorm(mu=0.05, sigma=2.0, N=300.),
    ...                           bins=200, kappa=0.6)

    Monodisperse sodium chloride:

    >>> aerosol2 = AerosolSpecies('NaCl', {'r_drys': [0.25], 'Nis': [1000.0]},
    ...                           kappa=0.2)
    """

    def __init__(
        self,
        species: str,
        distribution: BaseDistribution | dict[str, Any],
        kappa: float,
        rho: float | None = None,
        mw: float | None = None,
        bins: int | NDArray[np.floating[Any]] | None = None,
        r_min: float | None = None,
        r_max: float | None = None,
    ) -> None:
        self.species = species
        self.kappa = kappa
        self.rho = rho
        self.mw = mw
        self.distribution = distribution

        if isinstance(distribution, dict):
            self.r_drys = np.array(distribution["r_drys"]) * 1e-6

            if len(self.r_drys) > 1:
                # Reconstruct bin edges: assume r_dry is the geometric mean of its
                # bin, and the right edge of the first bin is the geometric mean of
                # the two smallest radii.
                mid1 = np.sqrt(self.r_drys[0] * self.r_drys[1])
                left = (self.r_drys[0] ** 2.0) / mid1
                rs = [left, mid1]
                for r in self.r_drys[1:]:
                    rs.append(r**2.0 / rs[-1])
                self.rs = np.array(rs) * 1e6
            else:
                self.rs = None  # pyrefly: ignore[bad-assignment]  # monodisperse; no bin edges needed

            self.Nis = np.array(distribution["Nis"])

        elif isinstance(distribution, Lognorm):
            if bins is None:
                raise ValueError("Need to specify `bins` when passing a Lognorm distribution")
            self.rs, self.r_drys, self.Nis = _discretize_dist(
                distribution,
                bins,
                r_min,
                r_max,
                mu_lo=distribution.mu,
                sigma_lo=distribution.sigma,
                mu_hi=distribution.mu,
                sigma_hi=distribution.sigma,
            )

        elif isinstance(distribution, MultiModeLognorm):
            if bins is None:
                raise ValueError(
                    "Need to specify `bins` when passing a MultiModeLognorm distribution"
                )
            self.rs, self.r_drys, self.Nis = _discretize_dist(
                distribution,
                bins,
                r_min,
                r_max,
                mu_lo=distribution.mus[0],
                sigma_lo=distribution.sigmas[0],
                mu_hi=distribution.mus[-1],
                sigma_hi=distribution.sigmas[-1],
            )

        else:
            raise ValueError(f"Unsupported distribution type {type(distribution)!r}")

        # total_N in cm⁻³ (before unit conversion); Nis converted to m⁻³
        self.total_N = float(np.sum(self.Nis))
        self.Nis = self.Nis * 1e6
        self.nr = len(self.r_drys)

    def stats(self) -> dict[str, float]:
        """Compute useful statistics about this aerosol's size distribution.

        Returns
        -------
        dict
            Statistics from the underlying distribution (see
            :meth:`~pyrcel.distributions.Lognorm.stats`). If ``rho`` was
            supplied, mass-weighted quantities are added: ``total_mass``,
            ``mean_mass``, and ``specific_surface_area``.

        Raises
        ------
        ValueError
            If the stored ``distribution`` does not implement ``stats()``.
        """
        if not isinstance(self.distribution, BaseDistribution):
            raise ValueError(
                f"Cannot compute stats for distribution of type {type(self.distribution)!r}"
            )

        stats_dict = self.distribution.stats()

        if self.rho:
            stats_dict["total_mass"] = stats_dict["total_volume"] * self.rho
            stats_dict["mean_mass"] = stats_dict["mean_volume"] * self.rho
            stats_dict["specific_surface_area"] = (
                stats_dict["total_surface_area"] / stats_dict["total_mass"]
            )

        return stats_dict

    def __repr__(self) -> str:
        return (
            f"AerosolSpecies({self.species!r}, kappa={self.kappa}, "
            f"distribution={self.distribution!r})"
        )

    def __str__(self) -> str:
        return f"{self.species}: kappa={self.kappa:.3f}, N={self.total_N:.1f} cm⁻³, nr={self.nr}"
