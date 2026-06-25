"""JAX-native Abdul-Razzak & Ghan (2000) activation parameterization.

Reference
---------
Abdul-Razzak, H., and S. J. Ghan (2000), A parameterization of aerosol
activation: 2. Multiple aerosol types, J. Geophys. Res., 105(D5), 6837-6844,
doi:10.1029/1999JD901161.

Ghan, S. J. et al (2011) Droplet Nucleation: Physically-based Parameterization
and Comparative Evaluation, J. Adv. Model. Earth Syst., 3,
doi:10.1029/2011MS000074.
"""

from __future__ import annotations

import jax

jax.config.update("jax_enable_x64", True)

import jax.numpy as jnp  # noqa: E402
from jax import Array  # noqa: E402
from jax.typing import ArrayLike  # noqa: E402

from .. import constants as c  # noqa: E402
from ..thermo import dv, dv_cont, es, ka_cont, sigma_w  # noqa: E402
from ._common import _kohler_crit_approx, _lognormal_act  # noqa: E402, F401


def arg2000(
    V: ArrayLike,
    T: ArrayLike,
    P: ArrayLike,
    mus: ArrayLike,
    sigmas: ArrayLike,
    Ns: ArrayLike,
    kappas: ArrayLike,
    accom: float = c.ac,
) -> tuple[Array, Array, Array]:
    """JAX-native Abdul-Razzak & Ghan (2000) activation parameterization.

    A faithful JAX re-implementation of [pyrcel.legacy.activation.arg2000][].
    All floating-point operations use `jax.numpy` so the computation is
    fully traceable and differentiable via `jax.grad`.

    The non-unity accommodation correction follows Ghan et al. (2011), eq. (40).

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
        Geometric standard deviations (> 1).
    Ns : array-like, shape (n_modes,)
        Total number concentrations (cm⁻³).
    kappas : array-like, shape (n_modes,)
        Hygroscopicity parameters.
    accom : float, optional
        Condensation accommodation coefficient (default `pyrcel.constants.ac`).

    Returns
    -------
    smax : float
        Maximum parcel supersaturation (decimal).
    N_acts : jnp.ndarray, shape (n_modes,)
        Activated number concentration per mode (cm⁻³).
    act_fracs : jnp.ndarray, shape (n_modes,)
        Activated number fraction per mode.

    See Also
    --------
    pyrcel.legacy.activation.arg2000 : NumPy reference implementation.
    """
    mus = jnp.asarray(mus, dtype=float)
    sigmas = jnp.asarray(sigmas, dtype=float)
    Ns = jnp.asarray(Ns, dtype=float)
    kappas = jnp.asarray(kappas, dtype=float)

    # Thermodynamic coefficients
    wv_sat = es(T - 273.15)
    alpha = (c.g * c.Mw * c.L) / (c.Cp * c.R * T**2) - (c.g * c.Ma) / (c.R * T)
    gamma = (c.R * T) / (wv_sat * c.Mw) + (c.Mw * c.L**2) / (c.Cp * c.Ma * T * P)

    # Reference condensation growth coefficient G (accom = 1)
    G_a0 = (c.rho_w * c.R * T) / (wv_sat * dv_cont(T, P) * c.Mw)
    G_b = (c.L * c.rho_w * (c.L * c.Mw / (c.R * T) - 1.0)) / (ka_cont(T) * T)
    G_0 = 1.0 / (G_a0 + G_b)

    # Kelvin coefficient (Köhler A)
    A_koh = (2.0 * sigma_w(T) * c.Mw) / (c.rho_w * c.R * T)

    n_modes = len(mus)
    Smis = []
    Sparts = []

    for i in range(n_modes):
        am = mus[i] * 1e-6  # μm → m
        sig = sigmas[i]
        N_si = Ns[i] * 1e6  # cm⁻³ → m⁻³
        kap = kappas[i]

        r_crit, Smi = _kohler_crit_approx(T, am, kap)

        # Scale G for non-unity accommodation (Ghan et al. 2011, eq. 40)
        if accom == 1.0:
            G = G_0
        else:
            G_a_ac = (c.rho_w * c.R * T) / (wv_sat * dv(T, r_crit, P, accom) * c.Mw)
            G_ac = 1.0 / (G_a_ac + G_b)
            G_a_ac1 = (c.rho_w * c.R * T) / (wv_sat * dv(T, r_crit, P, 1.0) * c.Mw)
            G_ac1 = 1.0 / (G_a_ac1 + G_b)
            G = G_0 * G_ac / G_ac1

        fi = 0.5 * jnp.exp(2.5 * jnp.log(sig) ** 2)
        gi = 1.0 + 0.25 * jnp.log(sig)

        zeta = (2.0 / 3.0) * A_koh * jnp.sqrt(alpha * V / G)
        etai = (alpha * V / G) ** 1.5 / (N_si * gamma * c.rho_w * 2.0 * jnp.pi)

        Spa = fi * (zeta / etai) ** 1.5
        Spb = gi * ((Smi**2 / (etai + 3.0 * zeta)) ** 0.75)
        Sparts.append((1.0 / Smi**2) * (Spa + Spb))
        Smis.append(Smi)

    smax = 1.0 / jnp.sqrt(jnp.sum(jnp.stack(Sparts)))

    N_acts = []
    act_fracs = []
    for i in range(n_modes):
        N_act, act_frac = _lognormal_act(smax, sigmas[i], Ns[i], Smis[i])
        N_acts.append(N_act)
        act_fracs.append(act_frac)

    return smax, jnp.stack(N_acts), jnp.stack(act_fracs)


class ARG2000:
    """Abdul-Razzak & Ghan (2000) activation scheme.

    A thin callable wrapper around `arg2000` that satisfies the
    [ActivationScheme][pyrcel.activation.ActivationScheme] interface.  Instantiate
    once; call repeatedly.

    Parameters
    ----------
    accom : float, optional
        Condensation accommodation coefficient (forwarded to every call).
    """

    def __init__(self, accom: float = c.ac) -> None:
        self.accom = accom

    def __call__(
        self,
        V: ArrayLike,
        T: ArrayLike,
        P: ArrayLike,
        mus: ArrayLike,
        sigmas: ArrayLike,
        Ns: ArrayLike,
        kappas: ArrayLike,
        accom: float | None = None,
    ) -> tuple[Array, Array, Array]:
        return arg2000(
            V, T, P, mus, sigmas, Ns, kappas, accom=self.accom if accom is None else accom
        )

    def __repr__(self) -> str:
        return f"ARG2000(accom={self.accom})"
