"""Unit tests for :mod:`pyrcel.console_report`."""

from __future__ import annotations

import numpy as np

from pyrcel.console_report import (
    PhaseTimer,
    equilibration_residual,
    print_summary,
    print_trajectory_table,
)


def test_equilibration_residual_zero_when_empty():
    assert equilibration_residual(280.0, 0.01, [], [], np.array([])) == 0.0


def test_phase_timer_logs_first_compile(capsys):
    calls = []

    def fn():
        calls.append(1)
        return np.array([1.0])

    timer = PhaseTimer()
    timer.run(("k",), "test phase", fn)
    timer.run(("k",), "test phase", fn)
    err = capsys.readouterr().err
    assert "Compiling test phase" in err
    assert err.count("finished in") == 2
    assert len(calls) == 2


def test_print_trajectory_table_subsamples(capsys):
    n = 100
    time = np.linspace(0, 10, n)
    x = np.zeros((n, 7))
    x[:, 0] = np.linspace(0, 1000, n)
    x[:, 2] = 280.0
    x[:, 4] = 1e-6
    x[:, 6] = 0.01
    print_trajectory_table(time, x, max_rows=10)
    out = capsys.readouterr().out
    assert "Trajectory (post-hoc sample)" in out
    assert out.count("\n") >= 10


def test_live_step_printer(capsys):
    from pyrcel.console_report import LiveStepPrinter

    printer = LiveStepPrinter()
    printer(1, 0.0, 0.0, 283.15, -2.0, 0.0, 0.0)
    printer.finish()
    out = capsys.readouterr().out
    assert "Integration loop" in out
    assert "283.15" in out
    assert "end of integration loop" in out


def test_print_summary(capsys):
    print_summary(
        {
            "S_max": 0.012,
            "t_smax": 42.0,
            "T_smax": 265.0,
            "z_smax": 1200.0,
            "per_species": [
                {
                    "species": "AS",
                    "eq_act_frac": 0.5,
                    "kn_act_frac": 0.4,
                    "N": 100.0,
                    "N_act": 50.0,
                }
            ],
            "total_act_frac": 0.5,
        }
    )
    out = capsys.readouterr().out
    assert "Simulation summary" in out
    assert "AS" in out
    assert "total activated fraction" in out
