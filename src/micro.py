"""
Microphysics constants and sub-routines for Nenes model
"""
import numpy as np

## Microphysics constants
g = 9.81 # Gravitational constant, m/s**2
Cp = 1004.0 # Specific heat of dry air at constant pressure, J/kg 
L = 2.5e6 # Latent heat of condensation, J/kg
rho_w = 1e3 # Density of water, kg/m**3
Rd = 287.0 # Gas constant for dry air, J/(kg K)
R = 8.314 # Universal gas constant, J/(mol K)
Mw = 18.0153/1e3 # Molecular weight of water, kg/mol
sigma_w = lambda T: 0.0761 - 1.55e-4*(T-273.15) # surface tension of water, J/m^2 given T in Kelvin

## NOT CORRECTING FOR NON-CONTINUUM EFFECTS
Dv = 0.3/1e4 # Diffusivity of water vapor in air, m^2/s
ka = lambda T: 419.*(5.69 + 0.017*(T-273.15))*1e-5 # thermal conductivty of air, W/(m K) given T in Kelvin

## Aerosol Constants
# Ammonium Sulfate
Ms = 0.13214 # Molecular weight, kg/mol
rho_s = 1.769*1e-3*1e6
rho_u = 1.769*1e-3*1e6
epsilon = 0.5 # mass fraction of soluble material in the dry particle
rho_p = rho_u/(1.-epsilon*(1.-(rho_u/rho_s))) # total wet particle density
nu = 3.0 # number of ions into which a solute dissolves



## AUXILIARY FUNCTIONS
def es(T):
    """Returns saturation vapor pressure (Pascal) at temperature T (Celsius)
    Formula 2.17 in Rogers&Yau"""
    return 611.2*np.exp(17.67*T/(T+243.5))

def Seq(T, r, r_dry, ns):
    '''Equilibrium supersaturation predicted by Kohler theory'''
    A = (2.*Mw*sigma_w(T))/(R*T*rho_w*r)
    B = (3.*ns*Mw*nu)/(4.*np.pi*rho_w*(r**3 - r_dry**3))
    return np.exp(A - B) - 1.0
