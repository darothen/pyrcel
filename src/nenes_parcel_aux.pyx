
cimport cython
cimport numpy as np
import numpy as np

from cython.parallel import prange
from libc.math cimport exp
from math import pi

## Import constants from nenes_parcel
cdef double Mw = 18.0153/1e3 # Molecular weight of water, kg/mol
cdef double R = 8.314 # Universal gas constant, J/(mol K)
cdef double rho_w = 1e3 # Density of water, kg/m**3
cdef double nu = 2.0
cdef Rd = 287.0 # Gas constant for dry air, J/(kg K)
cdef double g = 9.81 # Gravitational constant, m/s**2
cdef double Dv = 0.3/1e4 # Diffusivity of water vapor in air, m^2/s
cdef double L = 2.5e6 # Latent heat of condensation, J/kg
cdef double Cp = 1004.0 # Specific heat of dry air at constant pressure, J/kg 
cdef double PI = pi # Pi, constant

## AUXILIARY FUNCTIONS
cdef inline double sigma_w(double T) nogil:
    return 0.0761 - 1.55e-4*(T-273.15) # surface tension of water, J/m^2 given T in Kelvin

cdef inline double ka(double T): 
    return 419.*(5.69 + 0.017*(T-273.15))*1e-5 # thermal conductivty of air, W/(m K) given T in Kelvin

@cython.cdivision(True)
cdef inline double es(double T):
    """Returns saturation vapor pressure (Pascal) at temperature T (Celsius)
    Formula 2.17 in Rogers&Yau"""
    return 611.2*exp(17.67*T/(T+243.5))

@cython.cdivision(True)
cdef double Seq(double T, double r, double r_dry, double ns) nogil:
    '''Equilibrium supersaturation predicted by Kohler theory'''
    cdef double A, B, C
    A = (2.*Mw*sigma_w(T))/(R*T*rho_w*r)
    B = (3.*ns*Mw*nu)/(4.*PI*rho_w*(r**3 - r_dry**3))
    C = exp(A - B) - 1.0
    return C

## DERIVATIVE
def der(double t, np.ndarray[double, ndim=1] y, 
        int nr, np.ndarray[double, ndim=1] nss, np.ndarray[double, ndim=1] r_drys,
        np.ndarray[double, ndim=1] Nis, double V):
    return _der(t, y, nr, nss, r_drys, Nis, V)

@cython.cdivision(True)
@cython.boundscheck(False)
cdef np.ndarray[double, ndim=1] _der(double t, np.ndarray[double, ndim=1] y, 
                                     int nr, np.ndarray[double, ndim=1] nss, np.ndarray[double, ndim=1] r_drys,
                                     np.ndarray[double, ndim=1] Nis, double V):
    #print t
    
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
    cdef double dP_dt = (-g*P*V)/(Rd*T)
    
    # 2) dr_dt
    cdef double G_a = (rho_w*R*T)/(pv_sat*Dv*Mw)
    cdef double G_b = (L*rho_w*((L*Mw/(R*T))-1))/(ka(T)*T)
    cdef double G = 1./(G_a + G_b)
    
    cdef np.ndarray[double, ndim=1] drs_dt = np.empty(dtype="d", shape=(nr)) 
    cdef Py_ssize_t i
    cdef double r, r_dry, ns    

    for i in prange(nr, nogil=True):
        r = rs[i]
        r_dry = r_drys[i]
        ns = nss[i]
        drs_dt[i] = (G/r)*(S - Seq(T, r, r_dry, ns))
        
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
    S_a = P*dwv_dt/(0.622*pv_sat)
    S_b = 1.+S
    S_c = 0.622*L*dT_dt/(R*T**2) + g*V/(R*T)
    dS_dt = S_a - S_b*S_c
    
    cdef list x
    x = [dP_dt, dT_dt, dwv_dt, dwc_dt, dS_dt]
    x.extend(drs_dt)
    
    return np.array(x)
