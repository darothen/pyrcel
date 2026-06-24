"""Unit tests for AerosolSpecies (pyrcel.aerosol) and size distributions (pyrcel.distributions)."""

from __future__ import annotations

import numpy as np
import pytest
from hypothesis import given, settings
from hypothesis import strategies as st

from pyrcel.aerosol import AerosolSpecies
from pyrcel.distributions import Lognorm, MultiModeLognorm

# ---------------------------------------------------------------------------
# Helpers / fixtures
# ---------------------------------------------------------------------------

SULFATE = Lognorm(mu=0.05, sigma=2.0, N=300.0)
NACL = {"r_drys": [0.25], "Nis": [1000.0]}
TWO_MODE = MultiModeLognorm(mus=(0.02, 0.1), sigmas=(1.5, 2.0), Ns=(500.0, 200.0))


# ---------------------------------------------------------------------------
# AerosolSpecies — construction
# ---------------------------------------------------------------------------


def test_lognorm_construction_basic():
    aer = AerosolSpecies("(NH4)2SO4", SULFATE, kappa=0.6, bins=50)
    assert aer.species == "(NH4)2SO4"
    assert aer.kappa == 0.6
    assert aer.nr == 50
    assert aer.r_drys.shape == (50,)
    assert aer.Nis.shape == (50,)
    assert aer.total_N == pytest.approx(300.0, rel=1e-2)


def test_lognorm_construction_requires_bins():
    with pytest.raises(ValueError, match="bins"):
        AerosolSpecies("test", SULFATE, kappa=0.6)


def test_multimode_construction_basic():
    aer = AerosolSpecies("sea_salt", TWO_MODE, kappa=1.12, bins=100)
    assert aer.nr == 100
    assert aer.total_N == pytest.approx(700.0, rel=1e-2)


def test_multimode_construction_requires_bins():
    with pytest.raises(ValueError, match="bins"):
        AerosolSpecies("test", TWO_MODE, kappa=0.5)


def test_dict_construction_monodisperse():
    aer = AerosolSpecies("NaCl", NACL, kappa=0.2)
    assert aer.nr == 1
    assert aer.rs is None
    np.testing.assert_allclose(aer.r_drys[0], 0.25e-6)
    assert aer.total_N == pytest.approx(1000.0)


def test_dict_construction_multidisperse():
    d = {"r_drys": [0.1, 0.2, 0.4], "Nis": [100.0, 200.0, 50.0]}
    aer = AerosolSpecies("test", d, kappa=0.3)
    assert aer.nr == 3
    assert aer.rs is not None
    assert len(aer.rs) == 4  # nr + 1 bin edges
    assert aer.total_N == pytest.approx(350.0)


def test_unknown_distribution_raises():
    with pytest.raises(ValueError, match="Unsupported"):
        AerosolSpecies("bad", object(), kappa=0.5)


def test_nis_in_si_units():
    """Nis must be stored in m⁻³ (i.e. 1e6 × cm⁻³ input)."""
    aer = AerosolSpecies("test", NACL, kappa=0.5)
    # NACL has Nis=[1000] cm⁻³; stored value should be 1e9 m⁻³
    np.testing.assert_allclose(aer.Nis[0], 1000.0 * 1e6)


def test_bins_array_passthrough():
    """When bins is an array, it is used directly as bin edges."""
    edges = np.logspace(-3, 1, 21)  # 20 bins
    aer = AerosolSpecies("test", SULFATE, kappa=0.6, bins=edges)
    assert aer.nr == 20
    np.testing.assert_array_equal(aer.rs, edges)


def test_zero_kappa():
    aer = AerosolSpecies("insoluble", SULFATE, kappa=0.0, bins=10)
    assert aer.kappa == 0.0
    assert aer.nr == 10


# ---------------------------------------------------------------------------
# AerosolSpecies — __repr__ and __str__
# ---------------------------------------------------------------------------


def test_repr_contains_species_and_kappa():
    aer = AerosolSpecies("(NH4)2SO4", SULFATE, kappa=0.6, bins=50)
    r = repr(aer)
    assert "(NH4)2SO4" in r
    assert "kappa=0.6" in r
    assert "Lognorm" in r


def test_str_one_liner():
    aer = AerosolSpecies("NaCl", NACL, kappa=0.2)
    s = str(aer)
    assert "NaCl" in s
    assert "kappa" in s
    # Should be a single line
    assert "\n" not in s


# ---------------------------------------------------------------------------
# Lognorm — pdf / cdf / invcdf / moment
# ---------------------------------------------------------------------------


def test_pdf_positive():
    d = Lognorm(mu=0.05, sigma=1.5, N=500.0)
    x = np.logspace(-3, 0, 50)
    assert np.all(d.pdf(x) >= 0)


def test_cdf_monotone():
    d = Lognorm(mu=0.05, sigma=1.5, N=500.0)
    x = np.logspace(-3, 0, 50)
    cdf = d.cdf(x)
    assert np.all(np.diff(cdf) >= 0)


def test_cdf_bounds():
    d = Lognorm(mu=0.05, sigma=1.5, N=500.0)
    assert d.cdf(1e-8) == pytest.approx(0.0, abs=1e-6)
    assert d.cdf(1e3) == pytest.approx(500.0, rel=1e-4)


def test_invcdf_roundtrip():
    """invcdf(cdf(x)) ≈ x for interior values."""
    d = Lognorm(mu=0.05, sigma=1.5, N=1.0)
    x = np.array([0.01, 0.05, 0.1, 0.5])
    y = d.cdf(x)
    np.testing.assert_allclose(d.invcdf(y), x, rtol=1e-6)


def test_invcdf_roundtrip_large_N():
    """invcdf(cdf(x)) ≈ x when N >> 1 (cdf values not in [0,1])."""
    d = Lognorm(mu=0.05, sigma=1.5, N=300.0)
    x = np.array([0.02, 0.05, 0.1])
    np.testing.assert_allclose(d.invcdf(d.cdf(x)), x, rtol=1e-6)


def test_invcdf_validation():
    d = Lognorm(mu=0.05, sigma=1.5, N=100.0)
    with pytest.raises(ValueError):
        d.invcdf(-0.1)
    with pytest.raises(ValueError):
        d.invcdf(105.0)  # > N


def test_moment_zeroth():
    """0th moment equals total number concentration N."""
    d = Lognorm(mu=0.05, sigma=1.5, N=300.0)
    assert d.moment(0) == pytest.approx(300.0, rel=1e-12)


def test_moment_first():
    """1st moment equals N * mean_radius."""
    d = Lognorm(mu=0.05, sigma=1.5, N=1.0)
    # mean_radius = mu * exp(0.5 * ln(sigma)^2)
    expected = d.mu * np.exp(0.5 * np.log(d.sigma) ** 2)
    assert d.moment(1) == pytest.approx(expected, rel=1e-12)


def test_pdf_integrates_to_N():
    """Numerical integral of pdf over a wide range ≈ N."""
    d = Lognorm(mu=0.05, sigma=1.5, N=250.0)
    x = np.logspace(-4, 2, 5000)
    integral = np.trapezoid(d.pdf(x), x)
    assert integral == pytest.approx(250.0, rel=1e-3)


def test_lognorm_repr():
    d = Lognorm(mu=0.05, sigma=1.5, N=300.0)
    r = repr(d)
    assert "Lognorm" in r
    assert "mu" in r
    assert "sigma" in r


# ---------------------------------------------------------------------------
# Lognorm — property-based tests
# ---------------------------------------------------------------------------


@given(
    mu=st.floats(min_value=1e-4, max_value=10.0),
    sigma=st.floats(min_value=1.01, max_value=5.0),
    N=st.floats(min_value=1.0, max_value=1e5),
)
@settings(max_examples=200)
def test_pdf_nonnegative_hypothesis(mu, sigma, N):
    d = Lognorm(mu=mu, sigma=sigma, N=N)
    x = np.logspace(np.log10(mu) - 3, np.log10(mu) + 3, 20)
    assert np.all(d.pdf(x) >= 0)


@given(
    mu=st.floats(min_value=1e-4, max_value=10.0),
    sigma=st.floats(min_value=1.01, max_value=5.0),
    N=st.floats(min_value=1.0, max_value=1e5),
)
@settings(max_examples=200)
def test_cdf_monotone_hypothesis(mu, sigma, N):
    d = Lognorm(mu=mu, sigma=sigma, N=N)
    x = np.logspace(np.log10(mu) - 2, np.log10(mu) + 2, 20)
    assert np.all(np.diff(d.cdf(x)) >= 0)


@given(
    mu=st.floats(min_value=1e-4, max_value=10.0),
    sigma=st.floats(min_value=1.01, max_value=5.0),
    N=st.floats(min_value=1.0, max_value=1e5),
    k=st.integers(min_value=0, max_value=4),
)
@settings(max_examples=200)
def test_moment_positive_hypothesis(mu, sigma, N, k):
    d = Lognorm(mu=mu, sigma=sigma, N=N)
    assert d.moment(k) > 0


# ---------------------------------------------------------------------------
# AerosolSpecies — edge cases
# ---------------------------------------------------------------------------


def test_single_bin_lognorm():
    """bins=1 does not raise; shapes are correct.

    Note: total_N will be << N because a single bin captures very little of a
    broad Lognorm.  This is expected discretization behaviour, not a bug.
    """
    aer = AerosolSpecies("test", Lognorm(mu=0.05, sigma=1.5, N=100.0), kappa=0.6, bins=1)
    assert aer.nr == 1
    assert aer.r_drys.shape == (1,)
    assert aer.Nis.shape == (1,)
    assert aer.total_N > 0  # positive but not necessarily ≈ N for bins=1


def test_multi_mode_list_construction():
    """A list of two dict-spec aerosols round-trips through AerosolSpecies."""
    modes = [
        AerosolSpecies("acc", Lognorm(mu=0.05, sigma=2.0, N=300.0), kappa=0.6, bins=40),
        AerosolSpecies("coarse", Lognorm(mu=0.5, sigma=2.0, N=10.0), kappa=1.1, bins=10),
    ]
    assert modes[0].species == "acc"
    assert modes[1].species == "coarse"
    total_bins = sum(m.nr for m in modes)
    assert total_bins == 50


# ---------------------------------------------------------------------------
# MultiModeLognorm — property-based tests
# ---------------------------------------------------------------------------


@given(
    mu1=st.floats(min_value=1e-4, max_value=0.1),
    mu2=st.floats(min_value=0.1, max_value=10.0),
    sigma1=st.floats(min_value=1.01, max_value=3.0),
    sigma2=st.floats(min_value=1.01, max_value=3.0),
    N1=st.floats(min_value=1.0, max_value=1e4),
    N2=st.floats(min_value=1.0, max_value=1e4),
)
@settings(max_examples=100)
def test_multimode_cdf_monotone_hypothesis(mu1, mu2, sigma1, sigma2, N1, N2):
    d = MultiModeLognorm(mus=(mu1, mu2), sigmas=(sigma1, sigma2), Ns=(N1, N2))
    x = np.logspace(-4, 2, 30)
    assert np.all(np.diff(d.cdf(x)) >= 0)


@given(
    mu1=st.floats(min_value=1e-4, max_value=0.1),
    mu2=st.floats(min_value=0.1, max_value=10.0),
    sigma1=st.floats(min_value=1.01, max_value=3.0),
    sigma2=st.floats(min_value=1.01, max_value=3.0),
    N1=st.floats(min_value=1.0, max_value=1e4),
    N2=st.floats(min_value=1.0, max_value=1e4),
)
@settings(max_examples=100)
def test_multimode_total_N_is_sum(mu1, mu2, sigma1, sigma2, N1, N2):
    d = MultiModeLognorm(mus=(mu1, mu2), sigmas=(sigma1, sigma2), Ns=(N1, N2))
    # cdf at a very large radius integrates to the total number concentration
    assert d.cdf(1e6) == pytest.approx(N1 + N2, rel=1e-6)
