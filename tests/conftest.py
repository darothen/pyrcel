"""Shared pytest configuration and fixtures for the v2 validation harness.

Two responsibilities:

1. **Enable float64 in JAX before any array is created** (design doc §6.1). The
   parcel model is numerically hopeless in float32 (radii ~1e-8 m, the
   ``r**3 - r_dry**3`` cancellation in ``Seq``), so x64 must be on globally.
2. **Load the frozen oracle fixtures** baked from ``master`` by
   ``tools/oracle/generate_fixtures.py`` and expose them via pytest fixtures,
   parametrized over the shared scenario matrix in ``tests/scenarios.py``.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

TESTS_DIR = Path(__file__).resolve().parent
FIXTURE_DIR = TESTS_DIR / "fixtures" / "oracle"

# Make the dependency-light scenario module importable regardless of how pytest
# was invoked.
sys.path.insert(0, str(TESTS_DIR))

# Enable float64 in JAX as early as possible (before any jax array is created).
try:
    import jax

    jax.config.update("jax_enable_x64", True)
    _HAVE_JAX = True
except Exception:  # pragma: no cover - jax optional in some envs
    _HAVE_JAX = False

import scenarios  # noqa: E402

#: Number of bulk (non-radius) state variables; mirrors ``pyrcel.constants``.
N_STATE_VARS = 7


def have_jax() -> bool:
    return _HAVE_JAX


def load_fixture(name: str) -> dict:
    """Load a frozen oracle fixture by scenario name, or skip if absent."""
    path = FIXTURE_DIR / f"{name}.npz"
    if not path.exists():
        pytest.skip(f"oracle fixture {name!r} not generated (run generate_fixtures.py)")
    with np.load(path, allow_pickle=False) as data:
        return {k: data[k] for k in data.files}


@pytest.fixture(params=[s["name"] for s in scenarios.SCENARIOS])
def scenario_name(request) -> str:
    return request.param


@pytest.fixture
def oracle(scenario_name: str) -> dict:
    return load_fixture(scenario_name)
