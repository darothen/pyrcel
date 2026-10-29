"""Standalone regressions for this fix on unmodified upstream."""

import jax
import jax.numpy as jnp
import numpy as np
import pytest

from pyrcel import AerosolSpecies, ParcelModel
from pyrcel import constants as c
from pyrcel.activation import binned_activation
from pyrcel.integrator import find_smax


@pytest.mark.slow
@pytest.mark.parametrize("live", [False, True])
def test_event_peak_diagnostics_use_event_state(live):
    aerosol = AerosolSpecies("mono", {"r_drys": [0.05], "Nis": [100.0]}, 0.6)
    model = ParcelModel([aerosol], V=0.5, T0=283.15, S0=-0.01, P0=85000.0)
    t_peak, s_peak, y_peak, activated = find_smax(model.y0, model.args, 300.0)
    assert activated
    eq, kn, *_ = binned_activation(
        s_peak, y_peak[2], y_peak[7:], aerosol.r_drys, aerosol.Nis, aerosol.kappa
    )
    assert float(eq) == 1.0 and float(kn) == 1.0
    n_act = float(eq) * np.sum(aerosol.Nis)
    if len(model.args) == 6:
        # The sixth parameter fixes initial dry-air density. Compare the
        # concentration at the peak using the independent ideal-gas relation.
        rho_peak = y_peak[1] / (c.Rd * y_peak[2] * (1 + y_peak[3] / c.epsilon))
        n_act *= float(rho_peak) / float(model.args[5])

    model.run(
        300.0,
        output_dt=300.0,
        terminate=True,
        terminate_depth=100.0,
        live=live,
        live_chunk_dt=300.0,
    )
    summary = model.summary()
    assert len(model.time) == 2
    assert summary["t_smax"] == pytest.approx(float(t_peak))
    assert summary["z_smax"] == pytest.approx(float(y_peak[0]))
    assert summary["T_smax"] == pytest.approx(float(y_peak[2]))
    assert summary["per_species"][0]["kn_act_frac"] == pytest.approx(float(kn))
    assert summary["per_species"][0]["N_act"] == pytest.approx(n_act)


def test_interpolated_peak_state_matches_returned_peak_time(monkeypatch):
    import pyrcel.integrator as integrator

    times = np.array([0.0, 1.0, 2.0])
    states = np.zeros((3, 7))
    states[:, 0] = 10.0 * times
    states[:, 6] = 1.0 - (times - 0.4) ** 2

    def rhs(t, y, args):
        return jnp.zeros_like(y).at[0].set(10.0).at[6].set(-2.0 * (t - 0.4))

    monkeypatch.setattr(integrator, "parcel_ode_sys", rhs)
    maximum, time, state = integrator._peak_from_trajectory(times, states, (0, 0, 0, 0, 0))
    assert time == pytest.approx(0.4)
    assert maximum == pytest.approx(1.0)
    assert state[0] == pytest.approx(4.0)
    assert state[6] == pytest.approx(maximum)


def test_refined_peak_is_not_below_a_saved_endpoint(monkeypatch):
    import pyrcel.integrator as integrator

    times = np.array([0.0, 1.0, 2.0])
    states = np.zeros((3, 7))
    states[:, 6] = [1.0, 0.5, 0.0]

    def rhs(t, y, args):
        return jnp.zeros_like(y).at[6].set(jnp.interp(t, times, jnp.array([-2.0, -3.0, -1.0])))

    monkeypatch.setattr(integrator, "parcel_ode_sys", rhs)
    maximum, time, state = integrator._peak_from_trajectory(times, states, (0, 0, 0, 0, 0))
    assert maximum == 1.0
    assert time == 0.0
    assert state[6] == maximum


def test_differentiable_peak_is_not_below_a_saved_endpoint(monkeypatch):
    from types import SimpleNamespace

    import diffrax as dfx

    import pyrcel.integrator as integrator

    times = jnp.array([0.0, 1.0, 2.0])
    states = jnp.zeros((3, 7)).at[:, 6].set(jnp.array([1.0, 0.5, 0.0]))

    def rhs(t, y, args):
        return jnp.zeros_like(y).at[6].set(jnp.interp(t, times, jnp.array([-2.0, -3.0, -1.0])))

    monkeypatch.setattr(integrator, "parcel_ode_sys", rhs)
    monkeypatch.setattr(
        integrator,
        "_solve",
        lambda y0, *args: SimpleNamespace(
            ys=states.at[0, 6].set(y0[6]), result=dfx.RESULTS.successful
        ),
    )
    maximum = integrator.max_supersaturation(states[0], (0, 0, 0, 0, 0), times)
    assert float(maximum) == 1.0
    gradient = jax.grad(
        lambda value: integrator.max_supersaturation(
            states[0].at[6].set(value), (0, 0, 0, 0, 0), times
        )
    )(jnp.array(1.0))
    assert float(gradient) == pytest.approx(1.0)
