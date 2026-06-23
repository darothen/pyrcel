"""End-to-end high-level model tests (design §7.3 full pipeline; Phase 6).

Unlike ``test_trajectory`` (which seeds the v2 integrator with a frozen ``y0``),
these run the *whole* v2 pipeline -- :class:`pyrcel.model.ParcelModel` does its own
``optimistix`` equilibration, the diffrax solve with event termination, and the activation
diagnostics -- and check the headline numbers against the CVode oracle:

* ``S_max`` relative tolerance ``<= 1e-3`` (event-localized, so grid-independent),
* per-species activated fraction absolute tolerance ``<= 1e-3``.

Plus output-format coverage, the summary dict, a time-varying ``V(t)`` run, and that the
``console=True`` path prints without error.
"""

from __future__ import annotations

import numpy as np
import pytest
import scenarios as scn

from pyrcel.model import ParcelModel
from pyrcel.updraft import InterpolatedUpdraft

pytestmark = pytest.mark.slow

S_MAX_RTOL = 1e-3
ACT_FRAC_ATOL = 1e-3

_MODEL_CACHE: dict[str, ParcelModel] = {}


def _model(scenario_name: str) -> ParcelModel:
    if scenario_name not in _MODEL_CACHE:
        sc = scn.get_scenario(scenario_name)
        ic = sc["initial"]
        run = sc["run"]
        m = ParcelModel(
            scn.build_aerosols(sc),
            V=ic["V"],
            T0=ic["T0"],
            S0=ic["S0"],
            P0=ic["P0"],
            accom=ic["accom"],
        )
        m.run(
            run["t_end"], run["output_dt"], terminate=True, terminate_depth=run["terminate_depth"]
        )
        _MODEL_CACHE[scenario_name] = m
    return _MODEL_CACHE[scenario_name]


def test_equilibration_matches_oracle_y0(oracle, scenario_name):
    m = _model(scenario_name)
    np.testing.assert_allclose(m.y0, oracle["y0"], rtol=1e-9, atol=1e-12)


def test_smax_matches_oracle(oracle, scenario_name):
    m = _model(scenario_name)
    rel = abs(m.summary()["S_max"] - float(oracle["smax"])) / abs(float(oracle["smax"]))
    assert rel <= S_MAX_RTOL, f"S_max rel {rel:.2e}"


def test_activated_fraction_matches_oracle(oracle, scenario_name):
    m = _model(scenario_name)
    per = m.summary()["per_species"]
    for i, p in enumerate(per):
        assert abs(p["eq_act_frac"] - float(oracle["act_eq"][i])) <= ACT_FRAC_ATOL, p["species"]
        assert abs(p["kn_act_frac"] - float(oracle["act_kn"][i])) <= ACT_FRAC_ATOL, p["species"]


def test_output_formats():
    sc = scn.get_scenario("simple_sulfate")
    ic, run = sc["initial"], sc["run"]
    m = ParcelModel(
        scn.build_aerosols(sc),
        V=ic["V"],
        T0=ic["T0"],
        S0=ic["S0"],
        P0=ic["P0"],
        accom=ic["accom"],
    )
    parcel, aerosol = m.run(run["t_end"], run["output_dt"], mode="full")
    assert list(parcel.columns) == ["z", "P", "T", "wv", "wc", "wi", "S"]
    assert set(aerosol) == {"sulfate"}
    assert aerosol["sulfate"].shape[1] == m._nr

    x, heights = m.run(run["t_end"], run["output_dt"], mode="arrays")
    assert x.shape[0] == heights.shape[0]
    assert x.shape[1] == 7 + m._nr

    smax = m.run(run["t_end"], run["output_dt"], mode="smax")
    assert isinstance(smax, float) and smax > 0

    with pytest.raises(ValueError):
        m.run(run["t_end"], mode="nonsense")


def test_summary_structure():
    m = _model("simple_sulfate")
    s = m.summary()
    assert set(s) >= {"S_max", "t_smax", "T_smax", "per_species", "total_act_frac"}
    assert s["S_max"] > 0 and s["t_smax"] > 0
    assert 0.0 <= s["total_act_frac"] <= 1.0


def test_time_varying_updraft_runs():
    sc = scn.get_scenario("simple_sulfate")
    ic = sc["initial"]
    Vt = InterpolatedUpdraft(ts=[0.0, 40.0, 200.0], vs=[0.4, 1.2, 1.2])
    m = ParcelModel(
        scn.build_aerosols(sc),
        V=Vt,
        T0=ic["T0"],
        S0=ic["S0"],
        P0=ic["P0"],
        accom=ic["accom"],
    )
    smax = m.run(150.0, 1.0, terminate=True, terminate_depth=10.0, mode="smax")
    assert np.isfinite(smax) and smax > 0


def test_live_and_progress_mutually_exclusive():
    sc = scn.get_scenario("mono")
    ic, run = sc["initial"], sc["run"]
    m = ParcelModel(
        scn.build_aerosols(sc),
        V=ic["V"],
        T0=ic["T0"],
        S0=ic["S0"],
        P0=ic["P0"],
        accom=ic["accom"],
    )
    with pytest.raises(ValueError, match="mutually exclusive"):
        m.run(
            run["t_end"],
            run["output_dt"],
            live=True,
            progress=True,
        )


def test_live_output(capsys):
    sc = scn.get_scenario("mono")
    ic, run = sc["initial"], sc["run"]
    m = ParcelModel(
        scn.build_aerosols(sc),
        V=ic["V"],
        T0=ic["T0"],
        S0=ic["S0"],
        P0=ic["P0"],
        accom=ic["accom"],
        console=False,
    )
    m.run(
        run["t_end"],
        run["output_dt"],
        terminate=True,
        terminate_depth=run["terminate_depth"],
        live=True,
        live_chunk_dt=20.0,
    )
    out = capsys.readouterr().out
    assert "Integration loop" in out
    assert "end of integration loop" in out
    assert "|" in out  # z/T/S column divider


def test_console_output(capsys, caplog):
    import logging

    caplog.set_level(logging.INFO, logger="pyrcel")
    sc = scn.get_scenario("mono")
    ic, run = sc["initial"], sc["run"]
    m = ParcelModel(
        scn.build_aerosols(sc),
        V=ic["V"],
        T0=ic["T0"],
        S0=ic["S0"],
        P0=ic["P0"],
        accom=ic["accom"],
        console=True,
    )
    m.run(run["t_end"], run["output_dt"], terminate=True, terminate_depth=run["terminate_depth"])
    out = capsys.readouterr().out
    assert "Configuration" in out
    assert "Equilibrated initial state" in out
    assert "Integration plan" in out
    assert "Trajectory (post-hoc sample)" in out
    assert "Simulation summary" in out
    assert "S_max" in out
    assert "total activated fraction" in out
    assert "Compiling S_max event solve" in caplog.text
    assert "Termination:" in caplog.text


def test_backward_compat_alias():
    """ParcelModelJAX must remain importable and be identical to ParcelModel."""
    from pyrcel.model import ParcelModel, ParcelModelJAX

    assert ParcelModelJAX is ParcelModel
