"""Event-based S_max termination parity (design §4.5 option A, Phase 4).

The v2 model locates the supersaturation maximum with a ``dS/dt = 0`` event and stops
``terminate_depth`` metres past it, reproducing ``master``'s ``terminate=True``
behaviour. We assert:

* ``S_max`` parity with CVode (relative tolerance ``<= 1e-3``),
* ``t_smax`` parity -- the root-found ``t_smax`` lands within one output step of the
  (output-cadence-quantized) CVode ``t_smax``,
* the run terminates at ``z_smax + terminate_depth`` and stops well before ``t_end``
  (i.e. the event actually saved work).
"""

from __future__ import annotations

import jax.numpy as jnp
import numpy as np
import pytest
import scenarios as scn

from pyrcel.integrator_diffrax import find_smax, integrate_parcel_terminated

pytestmark = pytest.mark.slow

S_MAX_RTOL = 1e-3

_CACHE: dict[str, tuple] = {}


def _terminated(oracle, name: str):
    if name not in _CACHE:
        run = scn.get_scenario(name)["run"]
        y0 = jnp.asarray(oracle["y0"])
        args = (
            jnp.asarray(oracle["r_drys"]),
            jnp.asarray(oracle["Nis"]),
            jnp.asarray(oracle["kappas"]),
            float(oracle["accom"]),
            float(oracle["V"]),
        )
        ts, ys, info = integrate_parcel_terminated(
            y0,
            args,
            run["t_end"],
            run["output_dt"],
            terminate_depth=run["terminate_depth"],
        )
        _CACHE[name] = (ts, ys, info, run)
    return _CACHE[name]


def test_run_succeeds_and_finite(oracle, scenario_name):
    _, ys, info, _ = _terminated(oracle, scenario_name)
    assert info["success"]
    assert info["activated"]
    assert np.all(np.isfinite(ys))


def test_smax_value_matches_cvode(oracle, scenario_name):
    _, _, info, _ = _terminated(oracle, scenario_name)
    rel = abs(info["smax"] - float(oracle["smax"])) / abs(float(oracle["smax"]))
    assert rel <= S_MAX_RTOL, f"S_max rel {rel:.2e}"


def test_tsmax_matches_cvode(oracle, scenario_name):
    _, _, info, run = _terminated(oracle, scenario_name)
    # CVode t_smax is quantized to the output cadence; the root-found value must land
    # within one such step.
    assert abs(info["t_smax"] - float(oracle["t_smax"])) <= run["output_dt"]


def test_terminates_at_expected_depth(oracle, scenario_name):
    ts, ys, info, run = _terminated(oracle, scenario_name)
    V = float(oracle["V"])
    depth = run["terminate_depth"]
    # Final altitude should sit terminate_depth above z_smax (z = V t, so the residual
    # is V * |t_smax_v2 - t_smax_ref|, bounded by V * output_dt).
    z_final = float(ys[-1, 0])
    z_smax_ref = float(oracle["z_smax"])
    assert abs((z_final - depth) - z_smax_ref) <= V * run["output_dt"] + 1e-6
    # The event genuinely stopped the run early.
    assert info["t_cutoff"] < run["t_end"]
    # S_max is interior to the (terminated) trajectory.
    assert 0 < int(np.argmax(ys[:, 6])) < len(ys) - 1


def test_no_activation_returns_false(oracle, scenario_name):
    # A horizon far too short to reach the maximum reports activated=False
    # (and does not spuriously trigger the event).
    y0 = jnp.asarray(oracle["y0"])
    args = (
        jnp.asarray(oracle["r_drys"]),
        jnp.asarray(oracle["Nis"]),
        jnp.asarray(oracle["kappas"]),
        float(oracle["accom"]),
        float(oracle["V"]),
    )
    t_short = 0.5  # s; all scenarios reach S_max well after this
    _, _, _, activated = find_smax(y0, args, t_short)
    assert activated is False
