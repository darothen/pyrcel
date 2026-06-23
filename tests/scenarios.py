"""Canonical parcel-model scenarios shared by the oracle generator and the v2 tests.

This module is intentionally dependency-light: it contains only plain Python data
describing each scenario, plus a small ``build_aerosols`` helper that imports
``pyrcel`` lazily. That way it can be imported from either:

  * the *oracle* environment (``master`` checkout, numba + Assimulo/CVode), used to
    bake golden reference fixtures, and
  * the *v2* environment (JAX + diffrax), used to run the new model and compare.

Both backends must construct identical initial conditions, so the single source of
truth for "what configurations do we test?" lives here.
"""

from __future__ import annotations

# --- Scenario matrix -------------------------------------------------------------
#
# Each scenario is a dict:
#   name      : unique identifier (also the fixture filename stem)
#   initial   : parcel initial conditions + updraft + accommodation coefficient
#   aerosols  : list of aerosol-mode specs (lognormal or monodisperse)
#   run       : integration controls passed to ParcelModel.run
#
# Aerosol mode specs:
#   {"species", "kind": "lognorm", "mu", "sigma", "N", "kappa", "bins"}
#   {"species", "kind": "mono",    "r_drys": [...], "Nis": [...], "kappa"}
#
# Keep these fast: prefer modest bin counts and rely on `terminate` to stop early.

SCENARIOS: list[dict] = [
    {
        "name": "simple_sulfate",
        "initial": {"T0": 283.15, "S0": -0.02, "P0": 85000.0, "V": 1.0, "accom": 1.0},
        "aerosols": [
            {
                "species": "sulfate",
                "kind": "lognorm",
                "mu": 0.05,
                "sigma": 2.0,
                "N": 1000.0,
                "kappa": 0.54,
                "bins": 50,
            },
        ],
        "run": {
            "t_end": 200.0,
            "output_dt": 1.0,
            "solver_dt": 10.0,
            "terminate": True,
            "terminate_depth": 10.0,
            "max_steps": 2000,
        },
    },
    {
        "name": "two_mode",
        "initial": {"T0": 283.15, "S0": -0.05, "P0": 85000.0, "V": 0.44, "accom": 1.0},
        "aerosols": [
            {
                "species": "sulfate",
                "kind": "lognorm",
                "mu": 0.15,
                "sigma": 1.2,
                "N": 1000.0,
                "kappa": 0.54,
                "bins": 100,
            },
            {
                "species": "mos",
                "kind": "lognorm",
                "mu": 0.02,
                "sigma": 1.05,
                "N": 50000.0,
                "kappa": 0.12,
                "bins": 50,
            },
        ],
        "run": {
            "t_end": 400.0,
            "output_dt": 1.0,
            "solver_dt": 10.0,
            "terminate": True,
            "terminate_depth": 10.0,
            "max_steps": 2000,
        },
    },
    {
        "name": "fast_updraft",
        "initial": {"T0": 283.15, "S0": -0.02, "P0": 85000.0, "V": 5.0, "accom": 1.0},
        "aerosols": [
            {
                "species": "sulfate",
                "kind": "lognorm",
                "mu": 0.05,
                "sigma": 2.0,
                "N": 1000.0,
                "kappa": 0.54,
                "bins": 50,
            },
        ],
        "run": {
            "t_end": 60.0,
            "output_dt": 0.5,
            "solver_dt": 5.0,
            "terminate": True,
            "terminate_depth": 10.0,
            "max_steps": 2000,
        },
    },
    {
        "name": "slow_updraft",
        "initial": {"T0": 283.15, "S0": -0.02, "P0": 85000.0, "V": 0.2, "accom": 1.0},
        "aerosols": [
            {
                "species": "sulfate",
                "kind": "lognorm",
                "mu": 0.05,
                "sigma": 2.0,
                "N": 1000.0,
                "kappa": 0.54,
                "bins": 50,
            },
        ],
        "run": {
            "t_end": 600.0,
            "output_dt": 2.0,
            "solver_dt": 20.0,
            "terminate": True,
            "terminate_depth": 10.0,
            "max_steps": 2000,
        },
    },
    {
        "name": "low_accom",
        "initial": {"T0": 290.15, "S0": -0.01, "P0": 95000.0, "V": 1.0, "accom": 0.1},
        "aerosols": [
            {
                "species": "test",
                "kind": "lognorm",
                "mu": 0.025,
                "sigma": 1.82,
                "N": 1500.0,
                "kappa": 0.54,
                "bins": 100,
            },
        ],
        "run": {
            "t_end": 200.0,
            "output_dt": 1.0,
            "solver_dt": 10.0,
            "terminate": True,
            "terminate_depth": 50.0,
            "max_steps": 2000,
        },
    },
    {
        "name": "mono",
        "initial": {"T0": 283.15, "S0": -0.02, "P0": 85000.0, "V": 1.0, "accom": 1.0},
        "aerosols": [
            {
                "species": "NaCl",
                "kind": "mono",
                "r_drys": [0.1, 0.25, 0.5],
                "Nis": [500.0, 200.0, 50.0],
                "kappa": 1.2,
            },
        ],
        "run": {
            "t_end": 100.0,
            "output_dt": 1.0,
            "solver_dt": 10.0,
            "terminate": True,
            "terminate_depth": 10.0,
            "max_steps": 2000,
        },
    },
]


def get_scenario(name: str) -> dict:
    """Return the scenario dict with the given ``name``."""
    for scn in SCENARIOS:
        if scn["name"] == name:
            return scn
    raise KeyError(f"unknown scenario {name!r}")


def build_aerosols(scenario: dict):
    """Construct the list of :class:`pyrcel.AerosolSpecies` for ``scenario``.

    ``pyrcel`` is imported lazily so this module stays importable in environments
    where the model's heavy backend dependencies are unavailable.
    """
    import pyrcel as pm

    aerosols = []
    for spec in scenario["aerosols"]:
        kind = spec["kind"]
        if kind == "lognorm":
            dist = pm.Lognorm(mu=spec["mu"], sigma=spec["sigma"], N=spec["N"])
            aer = pm.AerosolSpecies(spec["species"], dist, kappa=spec["kappa"], bins=spec["bins"])
        elif kind == "mono":
            dist = {"r_drys": list(spec["r_drys"]), "Nis": list(spec["Nis"])}
            aer = pm.AerosolSpecies(spec["species"], dist, kappa=spec["kappa"])
        else:
            raise ValueError(f"unknown aerosol kind {kind!r}")
        aerosols.append(aer)
    return aerosols
