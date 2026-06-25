"""Abstract base class for activation parameterizations."""

from __future__ import annotations

from abc import ABC, abstractmethod
from typing import TYPE_CHECKING

from .. import constants as c

if TYPE_CHECKING:
    from jax import Array
    from jax.typing import ArrayLike


class ActivationScheme(ABC):
    """Interface for droplet activation parameterizations.

    All schemes accept lognormal aerosol parameters and updraft conditions and
    return the maximum parcel supersaturation plus per-mode activated fractions.
    JAX-native subclasses (e.g. [ARG2000][pyrcel.activation.ARG2000]) are fully differentiable via
    `jax.grad`.
    """

    @abstractmethod
    def __call__(
        self,
        V: ArrayLike,
        T: ArrayLike,
        P: ArrayLike,
        mus: ArrayLike,
        sigmas: ArrayLike,
        Ns: ArrayLike,
        kappas: ArrayLike,
        accom: float = c.ac,
    ) -> tuple[Array, Array, Array]:
        """Predict maximum supersaturation and activated fractions.

        Parameters
        ----------
        V : float
            Updraft speed (m/s).
        T : float
            Parcel temperature (K).
        P : float
            Parcel pressure (Pa).
        mus : array-like, shape (n_modes,)
            Median dry radii of each lognormal mode (μm).
        sigmas : array-like, shape (n_modes,)
            Geometric standard deviations (unitless, > 1).
        Ns : array-like, shape (n_modes,)
            Total number concentrations (cm⁻³).
        kappas : array-like, shape (n_modes,)
            Hygroscopicity parameters (unitless).
        accom : float, optional
            Condensation accommodation coefficient.

        Returns
        -------
        smax : float
            Maximum parcel supersaturation (decimal, e.g. 0.005 for 0.5%).
        N_acts : array, shape (n_modes,)
            Activated number concentration per mode (cm⁻³).
        act_fracs : array, shape (n_modes,)
            Activated number fraction per mode (0–1).
        """
