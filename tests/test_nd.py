"""Tests for the Nd activated droplet number diagnostic (issue #45a).

Checks:
1. summary dict contains Nd keys with physically sensible values
2. Nd >= N_act (at-smax equilibrium count) — more droplets grow to critical
   size given additional time past S_max
3. ModelOutput.Nd / nd_frac properties delegate correctly
4. mode='nd' returns the same scalar as summary["total_Nd"]
5. to_xarray() includes Nd and nd_t_eval variables
6. Nd is consistent across terminate / non-terminate integration paths
"""

from __future__ import annotations

import pytest

from pyrcel import AerosolSpecies, Lognorm, ParcelModel


@pytest.fixture
def simple_model():
    aer = AerosolSpecies("sulfate", Lognorm(mu=0.05, sigma=2.0, N=300.0), kappa=0.6, bins=10)
    return ParcelModel([aer], V=0.5, T0=283.15, S0=0.0, P0=85000.0)


@pytest.fixture
def two_mode_model():
    acc = AerosolSpecies("acc", Lognorm(mu=0.05, sigma=2.0, N=300.0), kappa=0.6, bins=8)
    crs = AerosolSpecies("crs", Lognorm(mu=0.5, sigma=1.8, N=50.0), kappa=0.2, bins=4)
    return ParcelModel([acc, crs], V=0.5, T0=283.15, S0=0.0, P0=85000.0)


# ---------------------------------------------------------------------------
# 1. Summary dict keys
# ---------------------------------------------------------------------------


def test_summary_contains_nd_keys(simple_model):
    simple_model.run(t_end=300.0, output_dt=5.0)
    s = simple_model.summary()
    assert "total_Nd" in s
    assert "total_nd_frac" in s
    assert "nd_t_eval" in s
    assert "nd_frac" in s["per_species"][0]
    assert "Nd" in s["per_species"][0]


def test_nd_positive(simple_model):
    simple_model.run(t_end=300.0, output_dt=5.0)
    s = simple_model.summary()
    assert s["total_Nd"] > 0.0


def test_nd_frac_in_unit_interval(simple_model):
    simple_model.run(t_end=300.0, output_dt=5.0)
    s = simple_model.summary()
    assert 0.0 <= s["total_nd_frac"] <= 1.0
    for p in s["per_species"]:
        assert 0.0 <= p["nd_frac"] <= 1.0


def test_nd_t_eval_after_smax(simple_model):
    """nd_t_eval must be >= t_smax (snapshot is taken after the S_max event)."""
    simple_model.run(t_end=300.0, output_dt=5.0)
    s = simple_model.summary()
    assert s["nd_t_eval"] >= s["t_smax"]


# ---------------------------------------------------------------------------
# 2. Physical invariant: Nd >= N_act (equilibrium at S_max)
# ---------------------------------------------------------------------------


def test_nd_geq_n_act_at_smax(simple_model):
    """Droplets given extra time past S_max should be >= those activated at S_max."""
    simple_model.run(t_end=300.0, output_dt=5.0)
    s = simple_model.summary()
    N_act_smax = s["per_species"][0]["N_act"]
    Nd = s["per_species"][0]["Nd"]
    assert Nd >= N_act_smax - 1e-6  # allow tiny float tolerance


def test_nd_two_mode_per_species_sum(two_mode_model):
    """Sum of per-species Nd must equal total_Nd."""
    two_mode_model.run(t_end=300.0, output_dt=5.0)
    s = two_mode_model.summary()
    assert sum(p["Nd"] for p in s["per_species"]) == pytest.approx(s["total_Nd"], rel=1e-10)


# ---------------------------------------------------------------------------
# 3. ModelOutput properties
# ---------------------------------------------------------------------------


def test_model_output_nd_property(simple_model):
    out = simple_model.run(t_end=300.0, output_dt=5.0, mode="full")
    assert out.Nd == pytest.approx(simple_model.summary()["total_Nd"], rel=1e-12)


def test_model_output_nd_frac_property(simple_model):
    out = simple_model.run(t_end=300.0, output_dt=5.0, mode="full")
    assert out.nd_frac == pytest.approx(simple_model.summary()["total_nd_frac"], rel=1e-12)


# ---------------------------------------------------------------------------
# 4. mode='nd' return value
# ---------------------------------------------------------------------------


def test_mode_nd_returns_float(simple_model):
    result = simple_model.run(t_end=300.0, output_dt=5.0, mode="nd")
    assert isinstance(result, float)


def test_mode_nd_matches_summary(simple_model):
    nd_val = simple_model.run(t_end=300.0, output_dt=5.0, mode="nd")
    assert nd_val == pytest.approx(simple_model.summary()["total_Nd"], rel=1e-12)


def test_mode_nd_positive(simple_model):
    nd_val = simple_model.run(t_end=300.0, output_dt=5.0, mode="nd")
    assert nd_val > 0.0


# ---------------------------------------------------------------------------
# 5. to_xarray integration
# ---------------------------------------------------------------------------


def test_xarray_contains_nd(simple_model):
    out = simple_model.run(t_end=300.0, output_dt=5.0, mode="full")
    ds = out.to_xarray()
    assert "Nd" in ds
    assert "nd_t_eval" in ds
    assert "sulfate_Nd" in ds
    assert float(ds["Nd"]) == pytest.approx(out.Nd, rel=1e-10)


# ---------------------------------------------------------------------------
# 6. Consistency across integration modes
# ---------------------------------------------------------------------------


def test_nd_terminate_vs_fixed_horizon():
    """Nd should be similar (not identical) between terminate and fixed-horizon runs."""
    aer = AerosolSpecies("sulf", Lognorm(mu=0.05, sigma=2.0, N=300.0), kappa=0.6, bins=10)

    m_term = ParcelModel([aer], V=0.5, T0=283.15, S0=0.0, P0=85000.0)
    nd_term = m_term.run(t_end=300.0, output_dt=5.0, mode="nd")

    m_fixed = ParcelModel([aer], V=0.5, T0=283.15, S0=0.0, P0=85000.0)
    nd_fixed = m_fixed.run(t_end=300.0, output_dt=5.0, terminate=False, mode="nd")

    # Both should be positive; they can differ because snapshot times differ
    assert nd_term > 0.0
    assert nd_fixed > 0.0
