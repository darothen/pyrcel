import numba as nb
import numpy as np
from numba.pycc import CC

import pyrcel.constants as c

## Define double DTYPE
DTYPE = np.float64

PI = 3.14159265358979323846264338328
N_STATE_VARS = c.N_STATE_VARS

# AOT/numba stuff
auxcc = CC("parcel_aux_numba")
auxcc.verbose = True


## Auxiliary, single-value calculations with GIL released for derivative
## calculations
@nb.njit()
@auxcc.export("sigma_w", "f8(f8)")
def sigma_w(T):
    """See :func:`pyrcel.thermo.sigma_w` for full documentation"""
    return 0.0761 - (1.55e-4) * (T - 273.15)


@nb.njit()
@auxcc.export("ka", "f8(f8, f8, f8)")
def ka(T, r, rho):
    """See :func:`pyrcel.thermo.ka` for full documentation"""
    ka_cont = 1e-3 * (4.39 + 0.071 * T)
    denom = 1.0 + (ka_cont / (c.at * r * rho * c.Cp)) * np.sqrt(
        (2 * PI * c.Ma) / (c.R * T)
    )
    return ka_cont / denom


@nb.njit()
@auxcc.export("dv", "f8(f8, f8, f8, f8)")
def dv(T, r, P, accom):
    """See :func:`pyrcel.thermo.dv` for full documentation"""
    P_atm = P * 1.01325e-5  # Pa -> atm
    dv_cont = 1e-4 * (0.211 / P_atm) * ((T / 273.0) ** 1.94)
    denom = 1.0 + (dv_cont / (accom * r)) * np.sqrt((2 * PI * c.Mw) / (c.R * T))
    return dv_cont / denom


@nb.njit()
@auxcc.export("es", "f8(f8)")
def es(T):
    """See :func:`pyrcel.thermo.es` for full documentation"""
    return 611.2 * np.exp(17.67 * T / (T + 243.5))


@nb.njit()
@auxcc.export("Seq", "f8(f8, f8, f8)")
def Seq(r, r_dry, T, kappa):
    """See :func:`pyrcel.thermo.Seq` for full documentation."""
    A = (2.0 * c.Mw * sigma_w(T)) / (c.R * T * c.rho_w * r)
    B = 1.0
    if kappa > 0.0:
        B = (r**3 - (r_dry**3)) / (r**3 - (r_dry**3) * (1.0 - kappa))
    return np.exp(A) * B - 1.0


## RHS Derivative callback function
@nb.njit(parallel=True)
@auxcc.export("parcel_ode_sys", "f8[:](f8[:], f8, i4, f8[:], f8[:], f8, f8[:], f8)")
def parcel_ode_sys(y, t, nr, r_drys, Nis, V, kappas, accom):
    """Calculates the instantaneous time-derivative of the parcel model system.

    Given a current state vector `y` of the parcel model, computes the tendency
    of each term including thermodynamic (pressure, temperature, etc) and aerosol
    terms. The basic aerosol properties used in the model must be passed along
    with the state vector (i.e. if being used as the callback function in an ODE
    solver).

    Parameters
    ----------
    y : array_like
        Current state of the parcel model system,
            * y[0] = altitude, m
            * y[1] = Pressure, Pa
            * y[2] = temperature, K
            * y[3] = water vapor mass mixing ratio, kg/kg
            * y[4] = cloud liquid water mass mixing ratio, kg/kg
            * y[5] = cloud ice water mass mixing ratio, kg/kg
            * y[6] = parcel supersaturation
            * y[7:] = aerosol bin sizes (radii), m
    t : float
        Current simulation time, in seconds.
    nr : Integer
        Number of aerosol radii being tracked.
    r_drys : array_like
        Array recording original aerosol dry radii, m.
    Nis : array_like
        Array recording aerosol number concentrations, 1/(m**3).
    V : float
        Updraft velocity, m/s.
    kappas : array_like
        Array recording aerosol hygroscopicities.
    accom : float, optional (default=:const:`constants.ac`)
        Condensation coefficient.

    Returns
    -------
    x : array_like
        Array of shape (``nr``+7, ) containing the evaluated parcel model
        instaneous derivative.

    Notes
    -----
    This function is implemented using numba; it does not need to be just-in-
    time compiled in order ot function correctly, but it is set up ahead of time
    so that the internal loop over each bin growth term is parallelized.

    """
    z = y[0]
    P = y[1]
    T = y[2]
    wv = y[3]
    wc = y[4]
    wi = y[5]
    S = y[6]
    rs = y[N_STATE_VARS:]

    T_c = T - 273.15  # convert temperature to Celsius
    pv_sat = es(T_c)  # saturation vapor pressure
    wv_sat = wv / (S + 1.0)  # saturation mixing ratio
    Tv = (1.0 + 0.61 * wv) * T
    e = (1.0 + S) * pv_sat  # water vapor pressure

    ## Compute air densities from current state
    rho_air = P / c.Rd / Tv
    #: TODO - port to parcel.py
    rho_air_dry = (P - e) / c.Rd / T

    ## Begin computing tendencies
    dP_dt = -1.0 * rho_air * c.g * V
    dwc_dt = 0.0
    # drs_dt = np.empty(shape=(nr), dtype=DTYPE)
    drs_dt = np.empty_like(rs)

    for i in nb.prange(nr):
        r = rs[i]
        r_dry = r_drys[i]
        kappa = kappas[i]
        Ni = Nis[i]

        ## Non-continuum diffusivity/thermal conductivity of air near
        ## near particle
        dv_r = dv(T, r, P, accom)
        ka_r = ka(T, r, rho_air)

        ## Condensation coefficient
        G_a = (c.rho_w * c.R * T) / (pv_sat * dv_r * c.Mw)
        G_b = (c.L * c.rho_w * ((c.L * c.Mw / (c.R * T)) - 1.0)) / (ka_r * T)
        G = 1.0 / (G_a + G_b)

        ## Difference between ambient and particle equilibrium supersaturation
        Seq_r = Seq(r, r_dry, T, kappa)
        delta_S = S - Seq_r

        ## Size and liquid water tendencies
        dr_dt = (G / r) * delta_S
        dwc_dt += (
            Ni * r * r * dr_dt
        )  # Contribution to liq. water tendency due to growth
        drs_dt[i] = dr_dt

    dwc_dt *= 4.0 * PI * c.rho_w / rho_air_dry  # Hydrated aerosol size -> water mass
    # use rho_air_dry for mixing ratio definition consistency
    # No freezing implemented yet
    dwi_dt = 0.0

    ## MASS BALANCE CONSTRAINT
    dwv_dt = -1.0 * (dwc_dt + dwi_dt)

    ## ADIABATIC COOLING
    dT_dt = -c.g * V / c.Cp - c.L * dwv_dt / c.Cp

    dz_dt = V

    """ Alternative methods for calculation supersaturation tendency
    # Used eq 12.28 from Pruppacher and Klett in stead of (9) from Nenes et al, 2001
    #cdef double S_a, S_b, S_c, dS_dt
    #cdef double S_b_old, S_c_old, dS_dt_old
    #S_a = (S+1.0)

    ## NENES (2001)
    #S_b_old = dT_dt*wv_sat*(17.67*243.5)/((243.5+(Tv-273.15))**2.)
    #S_c_old = (rho_air*g*V)*(wv_sat/P)*((0.622*L)/(Cp*Tv) - 1.0)
    #dS_dt_old = (1./wv_sat)*(dwv_dt - S_a*(S_b_old-S_c_old))

    ## PRUPPACHER (PK 1997)
    #S_b = dT_dt*0.622*L/(Rd*T**2.)
    #S_c = g*V/(Rd*T)
    #dS_dt = P*dwv_dt/(0.622*es(T-273.15)) - S_a*(S_b + S_c)

    ## SEINFELD (SP 1998)
    #S_b = L*Mw*dT_dt/(R*T**2.)
    #S_c = V*g*Ma/(R*T)
    #dS_dt = dwv_dt*(Ma*P)/(Mw*es(T-273.15)) - S_a*(S_b + S_c)
    """

    ## GHAN (2011)
    alpha = (c.g * c.Mw * c.L) / (c.Cp * c.R * (T**2))
    alpha -= (c.g * c.Ma) / (c.R * T)
    gamma = (P * c.Ma) / (c.Mw * pv_sat)
    gamma += (c.Mw * c.L * c.L) / (c.Cp * c.R * T * T)
    dS_dt = alpha * V - gamma * dwc_dt

    # x = np.empty(shape=(nr+N_STATE_VARS), dtype='d')
    x = np.empty_like(y)
    x[0] = dz_dt
    x[1] = dP_dt
    x[2] = dT_dt
    x[3] = dwv_dt
    x[4] = dwc_dt
    x[5] = dwi_dt
    x[6] = dS_dt
    x[N_STATE_VARS:] = drs_dt[:]

    return x
