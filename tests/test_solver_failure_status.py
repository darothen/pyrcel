"""Standalone regressions for this fix on unmodified upstream."""

import jax.numpy as jnp
import numpy as np
import pytest

from pyrcel import AerosolSpecies, ParcelModel
from pyrcel.ensemble import smax_nact_ensemble
from pyrcel.integrator import find_smax, max_supersaturation, nd_from_parcel
from pyrcel.util import ParcelModelError


def mono():
    a = AerosolSpecies("mono", {"r_drys": [0.05], "Nis": [100.0]}, 0.6)
    return ParcelModel([a], V=0.5, T0=283.15, S0=-0.01, P0=85000.0)


def test_failed_event_and_ensemble_status():
    m = mono()
    with pytest.raises(ParcelModelError, match="failed"):
        find_smax(m.y0, m.args, 300.0, max_steps=1)
    result = smax_nact_ensemble(m.y0, *m.args[:4], [0.5], 300.0, max_steps=1)
    assert not result["activated"][0]
    assert not result["success"][0]
    assert np.isnan(result["S_max"][0]) and np.isnan(result["N_act"][0])
    # A completed horizon without an event is a different, non-failing outcome.
    short = smax_nact_ensemble(m.y0, *m.args[:4], [0.5], 0.001)
    assert short["success"][0] and not short["activated"][0]


@pytest.mark.parametrize("diagnostic", [max_supersaturation, nd_from_parcel])
def test_failed_scalar_diagnostics_raise(diagnostic):
    m = mono()
    end = jnp.array([0.0, 1.0, 300.0]) if diagnostic is max_supersaturation else 300.0
    with pytest.raises(RuntimeError, match="integration failed"):
        diagnostic(m.y0, m.args, end, max_steps=1)
