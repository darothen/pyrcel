"""Standalone regressions for this fix on unmodified upstream."""

import numpy as np
import pytest
from scipy.optimize import brentq

from pyrcel import AerosolSpecies, ParcelModel
from pyrcel import constants as c
from pyrcel.activation._common import _kohler_crit_approx
from pyrcel.thermo import sigma_w


def test_nd_excludes_unactivated_tail():
    a = AerosolSpecies("mix", {"r_drys": [0.05, 0.5], "Nis": [100.0, 200.0]}, 0.6)
    m = ParcelModel([a], 0.5, 283.15, -0.01, 85000.0)
    rc, _ = _kohler_crit_approx(283.15, a.r_drys, 0.6)
    m.time = np.array([0.0, 1.0])
    m.x = np.tile(m.y0, (2, 1))
    m.x[-1, 6] = 0.01
    m.x[-1, 7:] = [2 * float(rc[0]), a.r_drys[1]]
    s = m._compute_summary((0.01, 1.0))
    assert s["total_Nd"] == pytest.approx(100e6)
    assert s["total_nd_frac"] == pytest.approx(1 / 3)

    # Both bins meet the equilibrium criterion; only the first is a droplet.
    assert s["per_species"][0]["N_act"] == pytest.approx(300e6)
    assert s["total_Nd"] < s["per_species"][0]["N_act"]


@pytest.mark.parametrize("kappa", [0.001, 0.004])
def test_subcritical_low_kappa_particles_are_not_droplets(kappa):
    aerosol = AerosolSpecies("low-kappa", {"r_drys": [0.1], "Nis": [100.0]}, kappa)
    model = ParcelModel([aerosol], 0.5, 283.15, -0.01, 85000.0)
    model.time = np.array([0.0, 1.0])
    model.x = np.tile(model.y0, (2, 1))
    # An independent scalar derivative root brackets the full-curve maximum.
    rd = 1e-7
    a = 2 * c.Mw * float(sigma_w(283.15)) / (c.R * 283.15 * c.rho_w * rd)
    x = brentq(
        lambda x: 3 * kappa * x**4 - a * (x**3 - 1) * (x**3 - 1 + kappa), 1.0, 1e4, xtol=1e-14
    )
    assert model.y0[7] < rd * x
    summary = model._compute_summary((-0.01, 0.0))
    assert summary["per_species"][0]["N_act"] == 0.0
    assert summary["total_Nd"] == 0.0


@pytest.mark.slow
def test_low_kappa_ascent_stays_below_activation():
    aerosol = AerosolSpecies("low-kappa", {"r_drys": [0.1], "Nis": [100.0]}, 0.001)
    model = ParcelModel([aerosol], 0.5, 283.15, -0.01, 85000.0)
    model.run(1.0, output_dt=0.5, terminate=False)
    summary = model.summary()
    assert summary["S_max"] < 0.0
    assert summary["per_species"][0]["N_act"] == 0.0
    assert summary["total_Nd"] == 0.0
