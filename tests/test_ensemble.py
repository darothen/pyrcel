"""Vectorized updraft-velocity ensemble (follow-up).

Validates the ``vmap``-ed ``S_max``/``N_act`` ensemble: sampling, batch shapes, physical
monotonicity (faster updraft -> higher peak supersaturation and more activated droplets),
and that a single-member ensemble reproduces the oracle ``S_max`` for ``simple_sulfate``.
"""

from __future__ import annotations

import numpy as np
import pytest
import scenarios as scn

from pyrcel.ensemble import (
    run_updraft_ensemble,
    sample_gaussian_updrafts,
    smax_nact_ensemble,
)
from pyrcel.equilibrate_jax import equilibrate_initial_state


def _simple_arrays():
    sc = scn.get_scenario("simple_sulfate")
    ic = sc["initial"]
    aers = scn.build_aerosols(sc)
    r_drys, kappas, Nis = [], [], []
    for aer in aers:
        r_drys.extend(aer.r_drys)
        kappas.extend([aer.kappa] * aer.nr)
        Nis.extend(aer.Nis)
    r_drys = np.asarray(r_drys)
    kappas = np.asarray(kappas)
    Nis = np.asarray(Nis)
    y0 = np.asarray(equilibrate_initial_state(ic["T0"], ic["S0"], ic["P0"], r_drys, kappas, Nis))
    return ic, y0, r_drys, Nis, kappas, aers


def test_sample_gaussian_updrafts():
    v = sample_gaussian_updrafts(0.5, 0.2, 1000, seed=1, v_min=0.05)
    assert v.shape == (1000,)
    assert np.all(v >= 0.05)
    # reproducible
    v2 = sample_gaussian_updrafts(0.5, 0.2, 1000, seed=1, v_min=0.05)
    np.testing.assert_array_equal(v, v2)


@pytest.mark.slow
def test_smax_increases_with_updraft():
    ic, y0, r_drys, Nis, kappas, _ = _simple_arrays()
    Vs = np.array([0.2, 0.5, 1.0, 2.0, 4.0])
    res = smax_nact_ensemble(y0, r_drys, Nis, kappas, ic["accom"], Vs, t_end=600.0)
    assert res["S_max"].shape == (5,)
    assert np.all(res["activated"])
    # strictly monotone in V (physical)
    assert np.all(np.diff(res["S_max"]) > 0), res["S_max"]
    assert np.all(np.diff(res["N_act"]) >= 0), res["N_act"]
    # activated number bounded by the total population
    assert np.all(res["N_act"] <= np.sum(Nis) + 1.0)


@pytest.mark.slow
def test_single_member_matches_oracle():
    ic, y0, r_drys, Nis, kappas, _ = _simple_arrays()
    d = scn.get_scenario("simple_sulfate")  # noqa: F841 (ensures scenario exists)
    res = smax_nact_ensemble(y0, r_drys, Nis, kappas, ic["accom"], np.array([1.0]), t_end=400.0)
    from conftest import load_fixture

    oracle = load_fixture("simple_sulfate")
    rel = abs(float(res["S_max"][0]) - float(oracle["smax"])) / float(oracle["smax"])
    assert rel <= 1e-3, rel


@pytest.mark.slow
def test_run_updraft_ensemble_endtoend():
    _, _, _, _, _, aers = _simple_arrays()
    res = run_updraft_ensemble(
        aers,
        T0=283.15,
        S0=-0.02,
        P0=85000.0,
        mean=0.5,
        std=0.2,
        n=16,
        seed=0,
        z_cap=300.0,
    )
    for key in ("S_max", "N_act", "T_smax", "activated", "V"):
        assert res[key].shape == (16,)
    assert np.all(np.isfinite(res["S_max"]))
    assert np.all(res["N_act"] >= 0)
    assert res["t_end"] > 0
