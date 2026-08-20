"""Shared utility functions for JAX-native activation parameterizations.

Extension (2026)
-----------------
``_kohler_crit_approx``, ``lognormal_activation``, ``binned_activation``, and
``multi_mode_activation`` now accept an optional surface-tension override.
It is named ``surf_tension`` (not ``sigma``/``sigmas``) *deliberately* --
``sigma``/``sigmas`` already means "geometric standard deviation of a
lognormal mode" throughout this module (see ``lognormal_activation``'s
``sigma`` parameter, and ``sigmas`` in ``_arg2000.py``/``_mbn2014.py``).
Reusing that name for surface tension would silently collide with an
unrelated, already-meaningful parameter. Default (``surf_tension=None``)
reproduces the original pure-water behavior exactly.

Note: ``_arg2000.py`` (ARG2000) and ``_mbn2014.py`` (MBN2014) have *not*
been updated with the same override -- they are independent activation
parameterizations not used by ``pyrcel.model.ParcelModel`` (which only
calls ``binned_activation`` from this file), and updating them safely
requires resolving the same naming collision there too. Flagged as a
known remaining gap, not fixed here.
"""

from __future__ import annotations

from collections.abc import Sequence

import jax

jax.config.update("jax_enable_x64", True)

import jax.numpy as jnp  # noqa: E402
from jax import Array  # noqa: E402
from jax.scipy.special import erfc  # noqa: E402
from jax.typing import ArrayLike  # noqa: E402

from .. import constants as c  # noqa: E402
from ..thermo import sigma_w  # noqa: E402


def _lognormal_act(
    smax: ArrayLike, sigma: ArrayLike, N: ArrayLike, sgi: ArrayLike
) -> tuple[Array, Array]:
    r"""Activated fraction of one lognormal mode.

    Parameters
    ----------
    smax : float
        Maximum parcel supersaturation (decimal).
    sigma : float
        Geometric standard deviation.
    N : float
        Total number concentration (any unit; returned N_act is in the same unit).
    sgi : float
        Modal critical supersaturation (decimal).

    Returns
    -------
    N_act : float
        Activated number concentration.
    act_frac : float
        Activated fraction (0–1).

    Notes
    -----
    Integrates the lognormal number distribution above the critical dry size
    corresponding to $S_\mathrm{max}$:

    $$N_\mathrm{act} = \frac{N}{2}\,\mathrm{erfc}(u_i), \qquad
      u_i = \frac{2\ln(s_{g,i}/S_\mathrm{max})}{3\sqrt{2}\ln\sigma}$$

    where $s_{g,i}$ is the modal critical supersaturation.
    """
    ui = 2.0 * jnp.log(sgi / smax) / (3.0 * jnp.sqrt(2.0) * jnp.log(sigma))
    N_act = 0.5 * N * erfc(ui)
    return N_act, N_act / N


def _kohler_crit_approx(
    T: ArrayLike,
    r_dry: ArrayLike,
    kappa: ArrayLike,
    surf_tension: ArrayLike | None = None,
) -> tuple[Array, Array]:
    r"""Approximate critical radius and supersaturation (JAX-traceable, kappa > 0).

    Vectorises over ``r_dry`` / ``kappa`` for scalar ``T``.

    Parameters
    ----------
    T : float
        Ambient temperature, K.
    r_dry : array or float
        Dry particle radius, m.
    kappa : array or float
        Hygroscopicity parameter.
    surf_tension : float, optional
        Surface tension override, J/m². Defaults to ``None`` (pure water,
        via `pyrcel.thermo.sigma_w`), matching the original implementation.
        Named ``surf_tension`` rather than ``sigma`` to avoid colliding with
        this module's existing use of ``sigma``/``sigmas`` for geometric
        standard deviation.

    Returns
    -------
    r_crit : Array
        Critical wet radius (m).
    s_crit : Array
        Critical supersaturation (decimal).

    Notes
    -----
    From the approximate κ-Köhler equation $S_\mathrm{eq} \approx A/r -
    \kappa r_d^3/r^3$, setting $dS_\mathrm{eq}/dr = 0$ gives:

    $$r_\mathrm{crit} = \sqrt{\frac{3\kappa r_d^3}{A}}, \qquad
      s_\mathrm{crit} = \sqrt{\frac{4A^3}{27\kappa r_d^3}}$$

    where $A = 2 M_w \sigma(T) / (\rho_w R T)$ is the Kelvin parameter,
    using $\sigma(T) = \sigma_w(T)$ by default or ``surf_tension`` otherwise.
    """
    sigma_eff = sigma_w(T) if surf_tension is None else surf_tension
    A = (2.0 * c.Mw * sigma_eff) / (c.rho_w * c.R * T)
    r_crit = jnp.sqrt((3.0 * kappa * r_dry**3) / A)
    s_crit = jnp.sqrt((4.0 * A**3) / (27.0 * kappa * r_dry**3))
    # When r_crit_approx < r_dry the approximate Köhler curve has no maximum
    # above the dry radius; the full curve also has no physical activation threshold.
    # Return inf so kinetic activation comparisons (rs >= r_crits) are False.
    r_crit = jnp.where(r_crit < r_dry, jnp.inf, r_crit)
    return r_crit, s_crit


def lognormal_activation(
    smax: ArrayLike,
    mu: ArrayLike,
    sigma: ArrayLike,
    N: ArrayLike,
    kappa: ArrayLike,
    T: ArrayLike | None = None,
    sgi: ArrayLike | None = None,
    surf_tension: ArrayLike | None = None,
) -> tuple[Array, Array]:
    r"""Activated number and fraction for one or more lognormal modes.

    At least one of ``T`` or ``sgi`` must be supplied.  If ``sgi`` is given
    it is used directly; otherwise the modal critical supersaturation is
    computed from ``T``, ``mu``, and ``kappa`` via the approximate κ-Köhler
    formula.

    All scalar arguments broadcast naturally, so the function handles both
    single-mode scalars and multi-mode arrays in one call.

    Parameters
    ----------
    smax : float
        Maximum parcel supersaturation (decimal).
    mu : float or array-like
        Geometric mean dry radius of the mode(s), m.
    sigma : float or array-like
        Geometric standard deviation(s). (Not surface tension -- see
        ``surf_tension`` below for that.)
    N : float or array-like
        Total number concentration(s) (cm⁻³ or any consistent unit).
    kappa : float or array-like
        Hygroscopicity parameter(s).
    T : float, optional
        Parcel temperature (K).  Required when ``sgi`` is not supplied.
    sgi : float or array-like, optional
        Pre-computed modal critical supersaturation(s).  If given, ``T``
        and ``mu`` and ``kappa`` are not used for the critical-point
        calculation.
    surf_tension : float, optional
        Surface tension override, J/m², forwarded to `_kohler_crit_approx`
        when ``sgi`` is not pre-supplied. Defaults to ``None`` (pure water).

    Returns
    -------
    N_act : float or array
        Activated number concentration (same unit as ``N``).
    act_frac : float or array
        Activated fraction in [0, 1].

    Notes
    -----
    Uses `_lognormal_act` with the modal critical supersaturation obtained
    from `_kohler_crit_approx`.

    """
    if sgi is None:
        if T is None:
            raise ValueError("Either T or sgi must be provided.")
        _, sgi = _kohler_crit_approx(T, mu, kappa, surf_tension)
    return _lognormal_act(smax, sigma, N, sgi)


def binned_activation(
    Smax: ArrayLike,
    T: ArrayLike,
    rs: ArrayLike,
    r_drys: ArrayLike,
    Nis: ArrayLike,
    kappa: ArrayLike,
    surf_tension: ArrayLike | None = None,
) -> tuple[Array, Array, Array, Array]:
    r"""Equilibrium and kinetic activation statistics for a binned aerosol mode.

    A JAX-native re-implementation of
    `pyrcel.legacy.activation.binned_activation` that operates on raw
    arrays (no `AerosolSpecies` object) and is `jax.jit`-compatible.
    Critical radii and supersaturations are computed with the approximate
    κ-Köhler formula (`_kohler_crit_approx`).

    Parameters
    ----------
    Smax : float
        Environmental maximum supersaturation (decimal).
    T : float
        Environmental temperature (K).
    rs : array-like, shape (nr,)
        Current wet radii of the aerosol/droplet population (m).
    r_drys : array-like, shape (nr,)
        Dry particle radii (m).
    Nis : array-like, shape (nr,)
        Number concentration per size bin (m⁻³).
    kappa : float
        Hygroscopicity parameter for this mode.
    surf_tension : float, optional
        Surface tension override for this mode, J/m². Defaults to ``None``
        (pure water, via `pyrcel.thermo.sigma_w`), matching the original
        implementation exactly. Set this to the same value used for this
        species elsewhere (e.g. `pyrcel.aerosol.AerosolSpecies.sigma`) so
        the activation diagnostics are consistent with the actual simulated
        physics for a surfactant-modified species.

    Returns
    -------
    eq_frac : float
        Equilibrium-activated fraction: bins with $S_\text{crit} \le S_\text{max}$.
    kn_frac : float
        Kinetic-activated fraction: bins that have grown past their critical
        wet radius.
    alpha : float
        Ratio $N_\text{kn} / N_\text{eq}$ (kinetic limitation factor).
    phi : float
        Fraction of kinetically active bins that are still below their
        critical wet radius (unactivated tail).

    Notes
    -----
    Follows the approach of Nenes et al. (2001): the equilibrium count includes
    all bins whose critical supersaturation is below the environmental value;
    the kinetic count starts at the smallest bin that has already grown past its
    critical wet radius and includes all larger bins.

    References
    ----------
    **[Nenes2001]** Nenes, A., Ghan, S., Abdul-Razzak, H., Chuang, P. Y., &
    Seinfeld, J. H. (2001). Kinetic limitations on cloud droplet formation
    and impact on cloud albedo. *Tellus B*, **53**(2), 133–149.
    https://doi.org/10.1034/j.1600-0889.2001.d01-12.x

    """
    rs = jnp.asarray(rs, dtype=float)
    r_drys = jnp.asarray(r_drys, dtype=float)
    Nis = jnp.asarray(Nis, dtype=float)

    N_tot = jnp.sum(Nis)
    r_crits, s_crits = _kohler_crit_approx(T, r_drys, kappa, surf_tension)

    # Equilibrium activation: bins whose critical supersaturation is below Smax.
    N_eq = jnp.sum(jnp.where(Smax >= s_crits, Nis, 0.0))
    eq_frac = N_eq / N_tot

    # Kinetic activation: among equilibrium-activated bins only, find the
    # smallest that has grown past its critical wet radius.  Gating on
    # s_crits <= Smax prevents the approximate Köhler formula from triggering
    # false positives for tiny particles near the r_crit ≈ r_dry transition.
    is_r_large = (rs >= r_crits) & (s_crits <= Smax)
    any_large = jnp.any(is_r_large)
    first_large = jnp.argmax(is_r_large)  # 0 when all False (gated by any_large)
    kn_mask = (jnp.arange(Nis.shape[0]) >= first_large) & any_large
    N_kn = jnp.sum(jnp.where(kn_mask, Nis, 0.0))
    kn_frac = N_kn / N_tot

    # Droplets inside the kinetic window that haven't reached critical size.
    N_unact = jnp.sum(jnp.where((rs < r_crits) & kn_mask, Nis, 0.0))
    phi = jnp.where(N_kn > 0.0, N_unact / N_kn, 1.0)
    alpha = jnp.where(N_eq > 0.0, N_kn / N_eq, 0.0)

    return eq_frac, kn_frac, alpha, phi


def multi_mode_activation(
    Smax: ArrayLike,
    T: ArrayLike,
    rss: Sequence[ArrayLike],
    r_dryss: Sequence[ArrayLike],
    Niss: Sequence[ArrayLike],
    kappas: ArrayLike,
    surf_tensions: Sequence[ArrayLike | None] | None = None,
) -> tuple[list[Array], list[Array]]:
    r"""Activation statistics for a multi-mode binned aerosol population.

    Calls [binned_activation][pyrcel.activation.binned_activation] for each
    mode and collects the equilibrium and kinetic activated fractions.

    Parameters
    ----------
    Smax : float
        Environmental maximum supersaturation (decimal).
    T : float
        Environmental temperature (K).
    rss : sequence of array-like
        Current wet radii for each mode, each of shape ``(nr_i,)`` (m).
    r_dryss : sequence of array-like
        Dry radii for each mode, each of shape ``(nr_i,)`` (m).
    Niss : sequence of array-like
        Number concentrations per bin for each mode, each shape ``(nr_i,)``
        (m⁻³).
    kappas : array-like, shape (n_modes,)
        Hygroscopicity parameter for each mode.
    surf_tensions : sequence, shape (n_modes,), optional
        Per-mode surface tension override, J/m², or ``None`` entries (or the
        whole argument as ``None``) for pure water. Defaults to ``None``
        (every mode uses pure water), matching the original implementation.

    Returns
    -------
    eq_fracs : list of float
        Equilibrium-activated fraction per mode.
    kn_fracs : list of float
        Kinetic-activated fraction per mode.

    See Also
    --------
    pyrcel.activation.binned_activation : Per-mode computation.
    """
    kappas = jnp.asarray(kappas, dtype=float)
    n_modes = len(kappas)
    if surf_tensions is None:
        surf_tensions = [None] * n_modes
    eq_fracs: list[Array] = []
    kn_fracs: list[Array] = []
    for rs, r_drys, Nis, kappa, st in zip(rss, r_dryss, Niss, kappas, surf_tensions):
        eq, kn, _, _ = binned_activation(Smax, T, rs, r_drys, Nis, kappa, st)
        eq_fracs.append(eq)
        kn_fracs.append(kn)
    return eq_fracs, kn_fracs
