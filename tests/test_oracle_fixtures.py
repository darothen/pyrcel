"""Smoke tests for the frozen oracle fixtures.

These do not exercise the (not-yet-written) v2 model. They guard the *integrity*
of the golden reference data so that later equivalence tests have a trustworthy
foundation: correct shapes, finiteness, internal consistency, and a manifest that
matches the scenario matrix. They also assert the global float64 invariant.
"""

from __future__ import annotations

import json

import numpy as np
import pytest
import scenarios
from conftest import FIXTURE_DIR, N_STATE_VARS, have_jax


def test_manifest_matches_scenarios():
    manifest = json.loads((FIXTURE_DIR / "manifest.json").read_text())
    fixture_names = {s["name"] for s in manifest["scenarios"]}
    expected = {s["name"] for s in scenarios.SCENARIOS}
    assert expected <= fixture_names, expected - fixture_names
    # Provenance must be recorded.
    assert manifest["git_commit"] != "unknown"
    assert "cvode" in manifest["solver"]["backend"].lower()
    assert manifest["versions"].get("assimulo", "").startswith("3.")


def test_shapes_consistent(oracle, scenario_name):
    nr = int(oracle["nr"])
    width = N_STATE_VARS + nr
    assert oracle["y0"].shape == (width,)
    assert oracle["r_drys"].shape == (nr,)
    assert oracle["kappas"].shape == (nr,)
    assert oracle["Nis"].shape == (nr,)
    assert oracle["rhs_Y"].ndim == 2 and oracle["rhs_Y"].shape[1] == width
    assert oracle["rhs_dYdt"].shape == oracle["rhs_Y"].shape
    assert oracle["traj_X"].ndim == 2 and oracle["traj_X"].shape[1] == width
    assert oracle["traj_t"].shape[0] == oracle["traj_X"].shape[0]
    assert int(oracle["n_physical"]) + int(oracle["n_random"]) == oracle["rhs_Y"].shape[0]


def test_all_finite(oracle):
    for key in (
        "y0",
        "r_drys",
        "kappas",
        "Nis",
        "rhs_Y",
        "rhs_dYdt",
        "traj_t",
        "traj_X",
    ):
        assert np.all(np.isfinite(oracle[key])), f"non-finite values in {key!r}"


def test_smax_consistent_with_trajectory(oracle):
    s_idx = N_STATE_VARS - 1  # S is the last bulk state var (index 6)
    smax_from_traj = float(oracle["traj_X"][:, s_idx].max())
    assert np.isclose(smax_from_traj, float(oracle["smax"]), rtol=0, atol=1e-12)


def test_radii_at_or_above_dry(oracle):
    """Every wet radius in the RHS sample set must satisfy r >= r_dry."""
    rs = oracle["rhs_Y"][:, N_STATE_VARS:]
    r_drys = oracle["r_drys"][None, :]
    # randomized states are clamped to r_dry * 1.0001; allow a hair of slack.
    assert np.all(rs >= r_drys * 0.999)


def test_initial_state_matches_first_trajectory_point(oracle):
    # The solver stores the initial condition as the first output row.
    np.testing.assert_allclose(oracle["traj_X"][0], oracle["y0"], rtol=1e-6, atol=1e-12)


def test_activation_fractions_in_range(oracle):
    for key in ("act_eq", "act_kn"):
        vals = oracle[key]
        assert np.all(vals >= -1e-12) and np.all(vals <= 1.0 + 1e-9), key


@pytest.mark.skipif(not have_jax(), reason="jax not installed in this environment")
def test_x64_enabled():
    import jax.numpy as jnp

    assert jnp.array(1.0).dtype == jnp.float64
    assert jnp.zeros(3).dtype == jnp.float64
