"""
.. module:: parcel
    :synopsis: Parcel model performance functions.

.. moduleauthor:: Daniel Rothenberg <darothen@mit.edu>

"""

from scipy.optimize import fminbound
import numpy as np

from numba import jit, float_, int32, int64

## Import constants from `micro` module
from micro import Mw, Ma, R, rho_w, Rd, g, Dv, ac, Ka, at, L, Cp
from micro import sigma_w, ka, dv, ka_T, Dv_T, es, Seq

## DERIVATIVE
@jit(float_[:](float_[:], float_, int32, float_[:], float_[:], float_, float_[:]))
def der(y, t, nr, r_drys, Nis, V, kappas):
    """Time-derivatives of variables tracked by the parcel model.

    :param y:
        Vector representing the current state of the parcel model system,
            * P - pressure (Pa)
            * T - temperature (K)
            * wv - water vapor mass mixing ratio (kg/kg)
            * wc - droplet liquid water mass mixing ratio (kg/kg)
            * S - parcel supersaturation
            * rs - droplet radii (m)
                -- the final *nr* elements of the state vector `y`
    :type y: np.ndarray

    :param t:
        Current time-step at which derivatives are being computed
    :type t: float

    :param nr:
        Total number of aerosol radii (across all species) being tracked within
        the parcel model and in the state vector `y`
    :type nr: int

    :param r_drys:
        Dry radii (m) of all the wetted aerosol being tracked, concatenated across all
        species into a single array of length `nr`. Should correspond to the order in
        which the wetted radii appear in `y`
    :type r_drys: np.ndarray

    :param Nis:
        Number concentrations (m**-3) of the wetted aerosols being tracked,
        concatenated across all species into a single array of length `nr`
    :type Nis: np.ndarray

    :param V:
        Updraft velocity (m)
    :type V: float

    :param kappas:
        Aerosol hygroscopicities, concatenated across all species into a single
        array of length `nr`
    :type kappas: np.ndarray

    :returns:
        derivatives of all variables in state vector `y`, in their original order

    """
    P = y[0]
    T = y[1]
    wv = y[2]
    wc = y[3]
    S = y[4]
    #P, T, wv, wc, S = y[:5]
    rs = np.array(y[5:])

    pv_sat = es(T-273.15) # saturation vapor pressure
    #wv_sat = wv/(S+1.) # saturation mixing ratio
    Tv = (1.+0.61*wv)*T # virtual temperature given parcel humidity
    rho_air = P/(Rd*Tv) # current air density accounting for humidity

    ## Calculate tendency terms
    # 1) Pressure
    dP_dt = (-g*P*V)/(Rd*Tv)

    # 2/3) Wet particle growth rates and droplet liquid water
    drs_dt = np.zeros(shape=(nr, ))
    dwc_dt = 0.
    for i in range(nr):
        r = rs[i]
        r_dry = r_drys[i]
        kappa = kappas[i]

        G_a = (rho_w*R*T)/(pv_sat*dv(T, r)*Mw)
        G_b = (L*rho_w*((L*Mw/(R*T))-1.))/(ka(T, rho_air, r)*T)
        G = 1./(G_a + G_b)

        ## Remove size-dependence from G
        #G_a = (rho_w*R*T)/(pv_sat*Dv*Mw)
        #G_b = (L*rho_w*((L*Mw/(R*T))-1.))/(Ka*T)

        delta_S = S - Seq(r, r_dry, T, kappa, 1.0)

        dr_dt = (G/r)*delta_S

        ## ---

        Ni = Nis[i]
        dwc_dt = dwc_dt + Ni*(r**2)*dr_dt
        drs_dt[i] = dr_dt
    dwc_dt = (4.*np.pi*rho_w/rho_air)*dwc_dt

    # 4) Water vapor content
    dwv_dt = -dwc_dt

    # 5) Temperature
    dT_dt = 0.
    dT_dt = -g*V/Cp - L*dwv_dt/Cp

    # 6) Supersaturation
    dS_dt = 0.
    ## NENES (2001) eq 9.
    #S_a = (S+1.0)
    #S_b_old = dT_dt*wv_sat*(17.67*243.5)/((243.5+(Tv-273.15))**2.)
    #S_c_old = (rho_air*g*V)*(wv_sat/P)*((0.622*L)/(Cp*Tv) - 1.0)
    #dS_dt_old = (1./wv_sat)*(dwv_dt - S_a*(S_b_old-S_c_old))

    ## PRUPPACHER (PK 1997) eq 12.28
    #S_a = (S+1.0)
    #S_b = dT_dt*0.622*L/(Rd*T**2.)
    #S_c = g*V/(Rd*T)
    #dS_dt = P*dwv_dt/(0.622*es(T-273.15)) - S_a*(S_b + S_c)

    ## SEINFELD (SP 1998)
    #S_a = (S+1.0)
    #S_b = L*Mw*dT_dt/(R*T**2.)
    #S_c = V*g*Ma/(R*T)
    #dS_dt = dwv_dt*(Ma*P)/(Mw*es(T-273.15)) - S_a*(S_b + S_c)

    ## GHAN (2011) - prefer to use this!
    alpha = (g*Mw*L)/(Cp*R*(T**2)) - (g*Ma)/(R*T)
    gamma = (P*Ma)/(Mw*es(T-273.15)) + (Mw*L*L)/(Cp*R*T*T)
    dS_dt = alpha*V - gamma*dwc_dt

    ## Repackage tendencies for feedback to numerical solver
    x = np.zeros(shape=(nr+5, ))
    x[0] = dP_dt
    x[1] = dT_dt
    x[2] = dwv_dt
    x[3] = dwc_dt
    x[4] = dS_dt
    x[5:] = drs_dt[:]

    ## Kill off unused variables to get rid of numba warnings
    extra = 0.*t*wc
    if extra > 1e6:
        print "used"

    return x
    #return 1.0
