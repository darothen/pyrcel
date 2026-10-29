"""Standalone regressions for this fix on unmodified upstream."""

import jax
import jax.numpy as jnp
import numpy as np
import pytest
from scipy.special import erfc

from pyrcel import constants as c
from pyrcel.activation import binned_activation
from pyrcel.activation._common import _kohler_crit_approx, _kohler_crit_exact, lognormal_activation
from pyrcel.equilibrate import kohler_crit
from pyrcel.thermo import Seq, sigma_w


def test_low_kappa_has_a_finite_physical_peak():
    rd, k, T = 1e-7, 0.001, 283.15
    rc, sc = _kohler_crit_exact(T, rd, k)
    assert float(rc) > rd and np.isfinite(float(rc))
    assert float(rc) == pytest.approx(float(kohler_crit(T, rd, k)), rel=1e-10)
    assert float(sc) == pytest.approx(float(Seq(rc, rd, T, k)), rel=1e-10)
    eq, kn, *_ = binned_activation(0.02, T, [2 * rc], [rd], [1e8], k)
    assert float(eq) == 1 and float(kn) == 1


def test_exact_fallback_is_differentiable_and_batched():
    def radius(k):
        return _kohler_crit_exact(283.15, 1e-7, k)[0]

    g = float(jax.jit(jax.grad(radius))(jnp.array(0.001)))
    h = 1e-7
    fd = float((radius(0.001 + h) - radius(0.001 - h)) / (2 * h))
    assert np.isfinite(g) and g == pytest.approx(fd, rel=1e-5)
    rc, _ = jax.jit(_kohler_crit_exact)(283.15, jnp.array([1e-7, 1e-5]), 0.001)
    assert np.isfinite(np.asarray(rc)).all()


def test_critical_point_has_no_approximation_boundary_jump():
    T, rd = 283.15, 1e-7
    A = 2 * c.Mw * float(sigma_w(T)) / (c.R * T * c.rho_w)
    boundary = A / (3 * rd)
    kappas = jnp.array([boundary * (1 - 1e-6), boundary * (1 + 1e-6)])
    radius, supersat = _kohler_crit_exact(T, rd, kappas)
    assert float(radius[1] / radius[0]) == pytest.approx(1.0, rel=2e-6)
    assert float(supersat[1] / supersat[0]) == pytest.approx(1.0, rel=2e-6)


def test_critical_points_are_maxima_of_full_water_activity_curve():
    temps, dry, k = np.meshgrid(
        [260.0, 283.15, 310.0],
        [5e-9, 1e-7, 1e-6],
        [1e-4, 0.001, 0.0038, 0.2, 0.6, 1.4],
        indexing="ij",
    )
    T, rd, kap = map(jnp.asarray, [temps.ravel(), dry.ravel(), k.ravel()])
    rc, sc = _kohler_crit_exact(T, rd, kap)
    deriv = jax.vmap(jax.grad(Seq, argnums=0))(rc, rd, T, kap)
    np.testing.assert_allclose(np.asarray(rc * deriv), 0.0, atol=1e-10)
    assert np.all(np.asarray(rc) > np.asarray(rd))
    assert np.all(np.asarray(sc) >= np.asarray(Seq(rc * 1.001, rd, T, kap)))
    assert np.all(np.asarray(sc) >= np.asarray(Seq(rc * 0.999, rd, T, kap)))


def test_lognormal_analytic_integral_keeps_its_power_law():
    # This covers the old hybrid's switch: s_crit in an analytic lognormal
    # inversion must retain r_dry**(-3/2), even where that approximation is poor.
    T, rd, kappa, sigma, smax = 283.15, 1e-7, 0.001, 1.8, 0.01
    A = 2 * c.Mw * float(sigma_w(T)) / (c.R * T * c.rho_w)
    cutoff = (4 * A**3 / (27 * kappa * smax**2)) ** (1 / 3)
    expected = 0.5 * erfc(np.log(cutoff / rd) / (np.sqrt(2) * np.log(sigma)))
    _, fraction = lognormal_activation(smax, rd, sigma, 100.0, kappa, T=T)
    assert float(fraction) == pytest.approx(expected, rel=1e-12)
    _, s1 = _kohler_crit_approx(T, rd, kappa)
    _, s2 = _kohler_crit_approx(T, 2 * rd, kappa)
    assert float(s2 / s1) == pytest.approx(2 ** (-1.5), rel=1e-12)


def test_kappa_water_activity_matches_solute_and_water_volumes():
    rd, r, kappa, T = 1e-7, 1.8e-7, 0.1, 283.15
    solute_volume = 4 * np.pi / 3 * rd**3
    water_volume = 4 * np.pi / 3 * (r**3 - rd**3)
    water_activity = 1 / (1 + kappa * solute_volume / water_volume)
    kelvin = np.exp(2 * c.Mw * float(sigma_w(T)) / (c.R * T * c.rho_w * r))
    assert float(Seq(r, rd, T, kappa)) + 1 == pytest.approx(water_activity * kelvin, rel=1e-13)


def test_ensemble_includes_exact_critical_threshold(monkeypatch):
    from types import SimpleNamespace

    import pyrcel.ensemble as ensemble

    y = jnp.array([0.0, 85000.0, 283.15, 0.009, 0.0, 0.0, 0.01, 1e-6])
    monkeypatch.setattr(
        ensemble,
        "_solve_to_smax",
        lambda *args: SimpleNamespace(ys=y[None, :], ts=jnp.array([10.0])),
    )
    monkeypatch.setattr(
        ensemble,
        "_kohler_crit_exact",
        lambda T, r, k: (jnp.ones_like(r), jnp.broadcast_to(y[6], r.shape)),
    )
    result = ensemble.smax_nact_ensemble(y, [5e-8], [1e8], [0.6], 1.0, [0.5], 300.0)
    assert result["activated"][0]
    assert result["N_act"][0] == pytest.approx(1e8)
