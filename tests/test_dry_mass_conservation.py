"""Standalone regressions for this fix on unmodified upstream."""

import jax
import jax.numpy as jnp
import numpy as np
import pytest

from pyrcel import AerosolSpecies, ParcelModel
from pyrcel import constants as c
from pyrcel.equilibrate import equilibrate_initial_state
from pyrcel.integrator import atol_vector, integrate_parcel, integrate_parcel_chunked
from pyrcel.parcel_aux import dry_air_density, parcel_args, parcel_ode_sys
from pyrcel.thermo import es


def mono():
    a = AerosolSpecies("mono", {"r_drys": [0.05], "Nis": [100.0]}, 0.6)
    return ParcelModel([a], V=0.5, T0=283.15, S0=-0.01, P0=85000.0)


def test_initial_water_definitions():
    T, P, RH, rd, N = 303.15, 80000.0, 0.5, 5e-8, 1e8
    y = np.asarray(equilibrate_initial_state(T, RH - 1, P, [rd], [0.6], [N]))
    e = RH * float(es(T - 273.15))
    assert y[3] == pytest.approx(c.epsilon * e / (P - e), rel=1e-13)
    liquid_mass_per_volume = 4 * np.pi / 3 * c.rho_w * N * (y[7] ** 3 - rd**3)
    assert y[4] == pytest.approx(liquid_mass_per_volume / ((P - e) / (c.Rd * T)), rel=1e-13)


def test_rhs_liquid_derivative_agrees_with_droplet_volume():
    m = mono()
    # Perturb pressure and humidity: this must not change number per dry mass.
    y = jnp.array(m.y0).at[1].set(70000.0).at[6].set(0.005)

    def liquid(state):
        return (
            4
            * jnp.pi
            / 3
            * c.rho_w
            * jnp.sum(m.args[1] / m.args[5] * (state[7:] ** 3 - m.args[0] ** 3))
        )

    dy = parcel_ode_sys(0.0, y, m.args)
    derivative = jax.jvp(liquid, (y,), (dy,))[1]
    assert float(derivative) == pytest.approx(float(dy[4]), rel=1e-12)


@pytest.mark.slow
def test_deep_ascent_conserves_droplet_water_and_chunk_reference():
    m = mono()
    out = m.run(2000.0, output_dt=100.0, terminate=False)
    # Check conservation at a precision tighter than the default radius atol;
    # dense interpolation of r^3 otherwise adds an O(1e-5) numerical residual.
    tight = integrate_parcel(m.y0, m.args, out.time, rtol=1e-10, atol=atol_vector(1) * 1e-3)
    ys = np.asarray(tight.ys)
    wc = (
        4
        * np.pi
        / 3
        * c.rho_w
        * np.sum(m.args[1] / m.args[5] * (ys[:, 7:] ** 3 - m.args[0] ** 3), axis=1)
    )
    np.testing.assert_allclose(ys[:, 4], wc, rtol=2e-6, atol=1e-10)
    ds = out.to_xarray()
    ni = ds["mono_number_concentration"].values[:, 0] * 1e6
    np.testing.assert_allclose(
        ni / np.asarray(dry_air_density(out.state)), m.args[1][0] / m.args[5], rtol=1e-12
    )
    args5 = m.args[:5]
    ts = np.linspace(0.0, 200.0, 5)
    fixed = integrate_parcel(m.y0, args5, ts)
    _, chunks, success = integrate_parcel_chunked(m.y0, args5, 200.0, 50.0, chunk_dt=50.0)
    assert success
    np.testing.assert_allclose(chunks, np.asarray(fixed.ys), rtol=2e-5, atol=1e-9)
    assert len(parcel_args(m.y0, args5)) == 6


def test_dry_mass_water_and_model_energy_tendencies():
    rd = jnp.array([5e-8, 2e-7])
    y = jnp.array([100.0, 80000.0, 281.0, 0.008, 0.001, 0.0, 0.003, 4e-7, 2e-6])
    args = (rd, jnp.array([1e8, 2e8]), jnp.array([0.6, 0.1]), 1.0, 0.5, 1.03)
    dy = parcel_ode_sys(30.0, y, args)
    assert float(dy[3] + dy[4] + dy[5]) == pytest.approx(0.0, abs=1e-20)
    # The invariant is for this model's constant-Cp/constant-L convention.
    dh = c.Cp * dy[2] + c.g * dy[0] + c.L * dy[3]
    assert float(dh) == pytest.approx(0.0, abs=1e-10)


def test_dry_air_density_uses_vapor_mixing_ratio():
    from pyrcel.parcel_aux import dry_air_density

    # S is prognostic and can disagree with the actual vapor mixing ratio.
    state = jnp.array([0.0, 70000.0, 275.0, 0.008, 0.001, 0.0, 0.02])
    expected = 70000.0 / (c.Rd * 275.0 * (1 + 0.008 / c.epsilon))
    assert float(dry_air_density(state)) == pytest.approx(expected, rel=1e-14)
    changed_S = state.at[6].set(-0.1)
    np.testing.assert_allclose(dry_air_density(changed_S), expected, rtol=1e-14)


@pytest.mark.parametrize("trim_initial", [False, True])
def test_exported_number_uses_vapor_and_original_reference(trim_initial):
    from pyrcel import AerosolSpecies
    from pyrcel.equilibrate import equilibrate_initial_state
    from pyrcel.model_output import ModelOutput

    aerosol = AerosolSpecies("mono", {"r_drys": [0.05], "Nis": [100.0]}, 0.6)
    T0, S0, P0 = 283.15, -0.01, 85000.0
    y0 = np.asarray(equilibrate_initial_state(T0, S0, P0, [5e-8], [0.6], [1e8]))
    end = y0.copy()
    end[1:4] = [70000.0, 275.0, 0.008]
    end[6] = 0.02
    states = np.stack([y0, end])
    times = np.array([0.0, 1000.0])
    summary = dict(
        S_max=0.02,
        t_smax=1000.0,
        total_Nd=0.0,
        nd_t_eval=1000.0,
        total_act_frac=0.0,
        per_species=[dict(species="mono", eq_act_frac=0.0, Nd=0.0)],
    )
    if trim_initial:
        states, times = states[1:], times[1:]
    out = ModelOutput(times, states, [aerosol], summary, 0.5, T0, S0, P0, 1.0)
    rho_initial = P0 / (c.Rd * T0 * (1 + y0[3] / c.epsilon))
    rho_end = 70000.0 / (c.Rd * 275.0 * (1 + 0.008 / c.epsilon))
    dataset = out.to_xarray()
    number = float(dataset["mono_number_concentration"].values[-1, 0])
    assert number == pytest.approx(100.0 * rho_end / rho_initial, rel=1e-13)
    assert float(dataset["rho"].values[-1]) == pytest.approx(rho_end * (1 + end[3]), rel=1e-13)


@pytest.mark.slow
def test_one_bin_ascent_matches_independent_scipy_reference():
    """Compare the coupled dry-mass trajectory with a NumPy/SciPy transcription."""
    from scipy.integrate import solve_ivp

    from pyrcel.thermo import dv_cont

    model = mono()
    rd, number, kappa, accom, speed, rho0 = model.args
    rd, number, kappa = map(np.asarray, (rd, number, kappa))
    times = np.array([0.0, 100.0, 400.0])
    # Share only the transport prefactor: pressure-unit corrections are tested
    # separately, and should not change what this conservation test asserts.
    reference_diffusivity = float(dv_cont(273.0, 101325.0))

    def reference_rhs(t, y):
        P, T, wv, S = y[1], y[2], y[3], y[6]
        radii = y[7:]
        tc = T - 273.15
        saturation = 611.2 * np.exp(17.67 * tc / (tc + 243.5))
        kelvin = 2 * c.Mw * (0.0761 - 1.55e-4 * tc) / (c.R * T * c.rho_w)
        water = (radii - rd) * (radii * radii + radii * rd + rd * rd)
        equilibrium = water / (water + kappa * rd**3) * np.exp(kelvin / radii) - 1
        density = P / (c.Rd * T * (1 + 0.61 * wv))
        diffusion = reference_diffusivity * (101325.0 / P) * (T / 273.0) ** 1.94
        diffusion /= 1 + diffusion / (accom * radii) * np.sqrt(2 * np.pi * c.Mw / (c.R * T))
        conductivity = 1e-3 * (4.39 + 0.071 * T)
        conductivity /= 1 + conductivity / (c.at * radii * density * c.Cp) * np.sqrt(
            2 * np.pi * c.Ma / (c.R * T)
        )
        growth = 1 / (
            c.rho_w * c.R * T / (saturation * diffusion * c.Mw)
            + c.L * c.rho_w * (c.L * c.Mw / (c.R * T) - 1) / (conductivity * T)
        )
        dr = growth / radii * (S - equilibrium)
        liquid = 4 * np.pi * c.rho_w * np.sum(number / rho0 * radii**2 * dr)
        alpha = c.g * c.Mw * c.L / (c.Cp * c.R * T**2) - c.g * c.Ma / (c.R * T)
        gamma = P * c.Ma / (c.Mw * saturation) + c.Mw * c.L**2 / (c.Cp * c.R * T**2)
        return np.r_[
            speed,
            -density * c.g * speed,
            -c.g * speed / c.Cp + c.L * liquid / c.Cp,
            -liquid,
            liquid,
            0.0,
            alpha * speed - gamma * liquid,
            dr,
        ]

    reference = solve_ivp(
        reference_rhs,
        (0.0, times[-1]),
        model.y0,
        method="BDF",
        t_eval=times,
        rtol=1e-10,
        atol=np.r_[np.full(3, 1e-8), np.full(4, 1e-13), 1e-16],
    )
    assert reference.success
    actual = np.asarray(
        integrate_parcel(
            model.y0,
            model.args,
            times,
            rtol=1e-10,
            atol=atol_vector(1) * 1e-3,
        ).ys
    )
    np.testing.assert_allclose(actual[:, 1:3], reference.y.T[:, 1:3], rtol=1e-8, atol=1e-7)
    np.testing.assert_allclose(actual[:, 3:7], reference.y.T[:, 3:7], rtol=2e-5, atol=1e-10)
    np.testing.assert_allclose(actual[:, 7:], reference.y.T[:, 7:], rtol=2e-5, atol=1e-12)


def test_summary_number_uses_each_snapshot_density():
    """A constructed expanded parcel can have lower Nd than its peak N_act."""
    model = mono()
    model.time = np.array([0.0, 1000.0])
    model.x = np.tile(model.y0, (2, 1))
    model.x[0, 6] = 0.02
    model.x[1, 1:4] = [70000.0, 275.0, 0.008]
    model.x[1, 6:] = [0.01, 1e-5]
    summary = model._compute_summary((0.02, 0.0))
    initial_density = model.y0[1] / (c.Rd * model.y0[2] * (1 + model.y0[3] / c.epsilon))
    final_density = 70000.0 / (c.Rd * 275.0 * (1 + 0.008 / c.epsilon))
    mode = summary["per_species"][0]
    assert mode["N_act"] == pytest.approx(100e6)
    assert mode["Nd"] == pytest.approx(100e6 * final_density / initial_density)
    assert mode["Nd"] < mode["N_act"]


@pytest.mark.parametrize("peak_form", ["state", "time", "trajectory"])
def test_summary_peak_density_matches_the_returned_state(peak_form, monkeypatch):
    import pyrcel.model as model_module

    model = mono()
    model.time = np.array([0.0, 1000.0])
    model.x = np.tile(model.y0, (2, 1))
    model.x[1, :4] = [500.0, 70000.0, 275.0, 0.008]
    model.x[1, 6:] = [0.01, 1e-5]
    peak_state = model.y0.copy()
    peak_state[:4] = [125.0, 81000.0, 282.0, 0.0088]
    peak_state[6:] = [0.02, 2e-6]
    if peak_form == "state":
        peak = (0.02, 250.0, peak_state)
    elif peak_form == "time":
        peak = (0.02, 250.0)
        peak_state = np.array([np.interp(250.0, model.time, column) for column in model.x.T])
    else:
        peak = None
        monkeypatch.setattr(
            model_module, "_peak_from_trajectory", lambda *args: (0.02, 250.0, peak_state)
        )
    summary = model._compute_summary(peak)
    rho_initial = model.y0[1] / (c.Rd * model.y0[2] * (1 + model.y0[3] / c.epsilon))
    rho_peak = peak_state[1] / (c.Rd * peak_state[2] * (1 + peak_state[3] / c.epsilon))
    assert summary["T_smax"] == pytest.approx(peak_state[2])
    assert summary["z_smax"] == pytest.approx(peak_state[0])
    assert summary["per_species"][0]["N_act"] == pytest.approx(100e6 * rho_peak / rho_initial)
