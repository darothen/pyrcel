"""Golden-trajectory regression: v2 diffrax vs frozen master/CVode (design §7.3).

Each scenario is integrated with the v2 :mod:`pyrcel.integrator_diffrax` solver,
**seeded from master's frozen ``y0``** and sampled at the frozen output times, then
compared to the CVode reference under the agreed acceptance tolerances:

* ``S_max`` relative tolerance ``<= 1e-3``
* activated-fraction absolute tolerance ``<= 1e-3``

(BDF and ESDIRK are different algorithms, so we assert *physical* agreement, not
bitwise equality -- §6.2.) Plus per-variable profile closeness and solver-agnostic
physics invariants.
"""

from __future__ import annotations

import jax.numpy as jnp
import numpy as np

import scenarios as scn
from pyrcel.activation import binned_activation
from pyrcel.integrator_diffrax import integrate_parcel_arrays

# Agreed acceptance tolerances (design §7, §9).
S_MAX_RTOL = 1e-3
ACT_FRAC_ATOL = 1e-3

# One solve per scenario, reused across the tests below.
_SOLUTION_CACHE: dict[str, tuple] = {}


def _solution(oracle, name: str):
    if name not in _SOLUTION_CACHE:
        y0 = jnp.asarray(oracle["y0"])
        ts = jnp.asarray(oracle["traj_t"])
        args = (
            jnp.asarray(oracle["r_drys"]),
            jnp.asarray(oracle["Nis"]),
            jnp.asarray(oracle["kappas"]),
            float(oracle["accom"]),
            float(oracle["V"]),
        )
        tt, ys, ok = integrate_parcel_arrays(y0, args, ts)
        _SOLUTION_CACHE[name] = (np.asarray(tt), np.asarray(ys), bool(ok))
    return _SOLUTION_CACHE[name]


def test_solver_succeeds_and_finite(oracle, scenario_name):
    _, ys, ok = _solution(oracle, scenario_name)
    assert ok, "diffrax solve did not report success"
    assert np.all(np.isfinite(ys)), "non-finite values in trajectory"


def test_smax_matches_cvode(oracle, scenario_name):
    _, ys, _ = _solution(oracle, scenario_name)
    smax_v2 = float(np.max(ys[:, 6]))
    smax_ref = float(oracle["smax"])
    rel = abs(smax_v2 - smax_ref) / abs(smax_ref)
    assert rel <= S_MAX_RTOL, f"S_max rel error {rel:.2e} > {S_MAX_RTOL}"


def test_activated_fraction_matches_cvode(oracle, scenario_name):
    _, ys, _ = _solution(oracle, scenario_name)
    si = int(np.argmax(ys[:, 6]))
    Smax = float(ys[si, 6])
    T_smax = float(ys[si, 2])
    rs = ys[si, 7:]

    aerosols = scn.build_aerosols(scn.get_scenario(scenario_name))
    offset = 0
    for i, aer in enumerate(aerosols):
        nr = aer.nr
        eq, kn, _, _ = binned_activation(Smax, T_smax, rs[offset : offset + nr], aer)
        offset += nr
        assert abs(eq - float(oracle["act_eq"][i])) <= ACT_FRAC_ATOL, aer.species
        assert abs(kn - float(oracle["act_kn"][i])) <= ACT_FRAC_ATOL, aer.species


def test_trajectory_profiles_match(oracle, scenario_name):
    _, ys, _ = _solution(oracle, scenario_name)
    X = np.asarray(oracle["traj_X"])
    g = np.all(np.isfinite(ys), axis=1)
    # Per-variable tolerances; observed agreement is 1-2 orders tighter than these.
    np.testing.assert_allclose(ys[g, 0], X[g, 0], atol=1e-4, rtol=0)  # z
    np.testing.assert_allclose(ys[g, 1], X[g, 1], rtol=1e-5, atol=1.0)  # P
    np.testing.assert_allclose(ys[g, 2], X[g, 2], atol=1e-3, rtol=0)  # T
    np.testing.assert_allclose(ys[g, 3], X[g, 3], atol=1e-8, rtol=0)  # wv
    np.testing.assert_allclose(ys[g, 4], X[g, 4], atol=1e-8, rtol=0)  # wc
    np.testing.assert_allclose(ys[g, 6], X[g, 6], atol=1e-5, rtol=0)  # S


def test_physical_invariants(oracle, scenario_name):
    _, ys, _ = _solution(oracle, scenario_name)
    # Monotonic ascent: z = ∫ V dt with V > 0.
    assert np.all(np.diff(ys[:, 0]) > 0), "altitude not monotonically increasing"
    # Total water wv + wc conserved (no ice, dwv = -dwc by construction).
    total_water = ys[:, 3] + ys[:, 4]
    assert np.max(np.abs(total_water - total_water[0])) < 1e-10, "water not conserved"
    # Supersaturation is single-peaked: the maximum is interior to the run.
    si = int(np.argmax(ys[:, 6]))
    assert 0 < si < len(ys) - 1, "S_max at trajectory boundary"
