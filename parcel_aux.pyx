# cython: embedsignature:True
# cython: profile=True
# (adds doc-strings accessible by Sphinx)
"""
.. module:: parcel
    :synopsis: Parcel model performance functions.

.. moduleauthor:: Daniel Rothenberg <darothen@mit.edu>

"""

from cython.parallel import prange
from scipy.optimize import fminbound
from libc.math cimport exp, sqrt
from math import pi

import numpy as np
cimport numpy as np
cimport cython

## Import constants from nenes_parcel
cdef double Mw = 18.0153/1e3 # Molecular weight of water, kg/mol
cdef double Ma = 28.9/1e3 # Molecular weight of dry air, kg/mol
cdef double R = 8.314 # Universal gas constant, J/(mol K)
cdef double rho_w = 1e3 # Density of water, kg/m**3
cdef double Rd = 287.0 # Gas constant for dry air, J/(kg K)
cdef double g = 9.81 # Gravitational constant, m/s**2
cdef double Dv = 3.e-5 # Diffusivity of water vapor in air, m^2/s
cdef double ac = 1.0 # condensation constant
cdef double Ka = 2.e-2 # Thermal conductivity of air, J/m/s/K
cdef double at = 0.96 # thermal accomodation coefficient
cdef double L = 2.5e6 # Latent heat of condensation, J/kg
cdef double Cp = 1004.0 # Specific heat of dry air at constant pressure, J/kg
cdef double PI = pi # Pi, constant

## AUXILIARY FUNCTIONS
cdef inline double sigma_w(double T) nogil:
    """Surface tension of water"""
    return 0.0761 - (1.55e-4)*(T-273.15) # surface tension of water, J/m^2 given T in Kelvin

@cython.cdivision(True)
cdef inline double ka(double T, double rho, double r) nogil:
    """Thermal conductivity of air, modified for non-continuum effects

    Revise with equation 17.71, Seinfeld and Pandis?"""
    cdef double denom, ka_cont
    ka_cont = 1e-3*(4.39 + 0.071*T)
    denom = 1.0 + (ka_cont/(at*r*rho*Cp))*sqrt((2*PI*Ma)/(R*T))
    return Ka/denom

@cython.cdivision(True)
cdef inline double dv(double T, double r) nogil:
    """Diffusivity of water vapor in air, modified for non-continuum effects

    Revise with equation 17.62, Seinfeld and Pandis?"""
    cdef double denom, dv_cont, P
    P = 1. # atm
    dv_cont = 1e-4*(0.211/P)*((T/273.)**1.94)
    denom = 1.0 + (Dv/(ac*r))*sqrt((2*PI*Mw)/(R*T))
    return Dv/denom

@cython.cdivision(True)
cdef inline double es(double T):
    """Returns saturation vapor pressure (Pascal) at temperature T (Celsius)
    Formula 2.17 in Rogers&Yau"""
    return 611.2*exp(17.67*T/(T+243.5))

@cython.cdivision(True)
cdef double Seq(double r, double r_dry, double T, double kappa) nogil:
    """See :func:`parcel_model.micro.Seq` for full documentation.
    """
    cdef double A = (2.*Mw*sigma_w(T))/(R*T*rho_w*r)
    cdef double B
    if kappa == 0.0:
        B = 1.0
    else:
        B = (r**3 - (r_dry**3))/(r**3 - (r_dry**3)*(1.-kappa))
    cdef double returnval = exp(A)*B - 1.0
    return returnval

@cython.cdivision(True)
cpdef double Seq_gil(double r, double r_dry, double T, double kappa):
    """Header of :func:`parcel_model.micro.Seq` ported to Cython and GIL-released.
    """
    return Seq(r, r_dry, T, kappa)

### GUESSING FUNCTION
@cython.boundscheck(False)
cpdef guesses(double T0, double S0, np.ndarray[double, ndim=1] r_drys,
              np.ndarray[double, ndim=1] kappas):
    """
    - DEPRECATED -
    As of 2/21/2013, a more accurate bisection method is used to compute the equilibrium
    wet radii for inactivate aerosols, and this routine is no longer necessary.

    Given a parcel's temperature and supersaturation as well as a size distribution
    of aerosols and their kappa parameters, computes seed guesses for determining
    the equilibrium wet radii of the inactivated aerosols or, when possible, the exact
    equilibrium solution.
    """
    cdef unsigned int nr = r_drys.shape[0]
    cdef np.ndarray[double, ndim=1] guesses = np.empty(dtype='d', shape=(nr))

    cdef np.ndarray[double, ndim=1] r_range, ss
    cdef double ri, rdi, ki
    cdef unsigned int i, j, idx_min
    for i in xrange(nr):
        rdi = r_drys[i]
        ki = kappas[i]
        r_range = np.arange(rdi+rdi/1000., rdi*10., rdi/1000.)
        ss = np.empty(dtype='d', shape=(len(r_range)))
        for j in prange(r_range.shape[0], nogil=True):
            ss[j] = Seq(r_range[j], rdi, T0, ki)
        idx_min = np.argmin(np.abs(ss - S0))
        guesses[i] = r_range[idx_min]

    return guesses

## DERIVATIVE
def der(np.ndarray[double, ndim=1] y, double t,
       int nr, np.ndarray[double, ndim=1] r_drys, np.ndarray[double, ndim=1] Nis,
       double V, np.ndarray[double, ndim=1] kappas):
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

    return _der(t, y, nr, r_drys, Nis, V, kappas)

@cython.cdivision(True)
@cython.boundscheck(False)
cdef np.ndarray[double, ndim=1] _der(double t, np.ndarray[double, ndim=1] y,
                                     int nr, np.ndarray[double, ndim=1] r_drys,
                                     np.ndarray[double, ndim=1] Nis, double V,
                                     np.ndarray[double, ndim=1] kappas):
    """Private function for computing the derivatives used to integrate
    the parcel model forward in time. See :func:`parcel_model.parcel_aux.der`
    for complete documentation.
    """
    cdef double P = y[0]
    cdef double T = y[1]
    cdef double wv = y[2]
    cdef double wc = y[3]
    cdef double S = y[4]
    cdef np.ndarray[double, ndim=1] rs = y[5:]

    cdef double pv_sat = es(T-273.15) # saturation vapor pressure
    cdef double wv_sat = wv/(S+1.) # saturation mixing ratio
    cdef double Tv = (1.+0.61*wv)*T
    cdef double rho_air = P/(Rd*Tv)

    # 1) dP_dt
    cdef double dP_dt = (-g*P*V)/(Rd*Tv)
    # FIX AT CONSTANT PRESSURE
    #dP_dt = 0.0

    # 2) dr_dt
    cdef double G_a, G_b, G
    cdef np.ndarray[double, ndim=1] drs_dt = np.empty(dtype="d", shape=(nr))
    cdef unsigned int i
    cdef double r, r_dry, delta_S, kappa

    for i in prange(nr, nogil=True):
        r = rs[i]
        r_dry = r_drys[i]
        kappa = kappas[i]

        ## Remove size-dependence from G
        G_a = (rho_w*R*T)/(pv_sat*dv(T, r)*Mw)
        #G_a = (rho_w*R*T)/(pv_sat*Dv*Mw)
        G_b = (L*rho_w*((L*Mw/(R*T))-1.))/(ka(T, rho_air, r)*T)
        #G_b = (L*rho_w*((L*Mw/(R*T))-1.))/(Ka*T)
        G = 1./(G_a + G_b)

        delta_S = S - Seq(r, r_dry, T, kappa)

        drs_dt[i] = (G/r)*delta_S

    # 3) dwc_dt
    cdef double dwc_dt = 0.0
    cdef double Ni, dr_dt

    for i in range(nr):
        Ni = Nis[i]
        r = rs[i]
        dr_dt = drs_dt[i]
        dwc_dt = dwc_dt + Ni*(r**2)*dr_dt
    dwc_dt = (4.*PI*rho_w/rho_air)*dwc_dt

    # 4) dwv_dt
    cdef double dwv_dt
    dwv_dt = -dwc_dt

    # 5) dT_dt
    cdef double dT_dt
    dT_dt = -g*V/Cp - L*dwv_dt/Cp

    # 6) dS_dt
    # Used eq 12.28 from Pruppacher and Klett in stead of (9) from Nenes et al, 2001
    cdef double S_a, S_b, S_c, dS_dt
    cdef double S_b_old, S_c_old, dS_dt_old
    S_a = (S+1.0)

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

    ## GHAN (2011)
    cdef double alpha, gamma
    alpha = (g*Mw*L)/(Cp*R*(T**2)) - (g*Ma)/(R*T)
    gamma = (P*Ma)/(Mw*es(T-273.15)) + (Mw*L*L)/(Cp*R*T*T)
    dS_dt = alpha*V - gamma*dwc_dt

    #print t, dS_dt, dS_dt_old

    cdef np.ndarray[double, ndim=1] x = np.empty(dtype='d', shape=(nr+5))
    x[0] = dP_dt
    x[1] = dT_dt
    x[2] = dwv_dt
    x[3] = dwc_dt
    x[4] = dS_dt
    x[5:] = drs_dt[:]

    return x
