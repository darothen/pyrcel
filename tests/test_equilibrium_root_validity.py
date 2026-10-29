"""Standalone regressions for this fix on unmodified upstream."""

import pytest
from scipy.optimize import brentq

from pyrcel import AerosolSpecies, ParcelModel
from pyrcel import constants as c
from pyrcel.equilibrate import equilibrate_radii, kohler_crit
from pyrcel.thermo import Seq, sigma_w


def test_zero_kappa_is_rejected_before_solving():
    with pytest.raises(ValueError, match="kappa > 0"):
        equilibrate_radii(283.15, -0.05, [5e-8], [0.0])
    bad = AerosolSpecies("insoluble", {"r_drys": [0.05], "Nis": [100.0]}, 0.0)
    with pytest.raises(ValueError, match="kappa > 0"):
        ParcelModel([bad], 0.5, 283.15, -0.05, 85000.0)


def test_equilibrium_failures_are_explicit():
    with pytest.raises(RuntimeError):
        equilibrate_radii(283.15, -0.05, [5e-8], [0.6], max_steps=1)
    with pytest.raises(RuntimeError, match="root is not contained"):
        equilibrate_radii(283.15, 0.5, [5e-8], [0.6])


@pytest.mark.parametrize(
    "T,rd,kappa,S",
    [(260.0, 5e-9, 0.001, -0.1), (283.15, 1e-7, 0.6, -0.01), (310.0, 1e-6, 1.4, -0.02)],
)
def test_scaled_roots_match_scipy(T, rd, kappa, S):
    a = 2 * c.Mw * float(sigma_w(T)) / (c.R * T * c.rho_w * rd)

    def derivative(x):
        delta = (x - 1) * (x * x + x + 1)
        return 3 * kappa * x**4 - a * delta * (delta + kappa)

    xc = brentq(derivative, 1.0, 1e4, xtol=1e-13)
    rc = float(kohler_crit(T, rd, kappa))
    assert rc / rd == pytest.approx(xc, rel=1e-11)
    xr = brentq(lambda x: float(Seq(x * rd, rd, T, kappa)) - S, 1.0, xc, xtol=1e-13)
    r = float(equilibrate_radii(T, S, [rd], [kappa])[0])
    assert r / rd == pytest.approx(xr, rel=1e-11)


@pytest.mark.parametrize("kappa", [-1.0, float("nan"), float("inf")])
def test_nonfinite_and_negative_kappa_are_rejected(kappa):
    with pytest.raises(ValueError, match="kappa > 0"):
        equilibrate_radii(283.15, -0.01, [5e-8], [kappa])
