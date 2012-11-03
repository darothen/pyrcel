"""
Microphysics constants and sub-routines for Nenes model
"""
import numpy as np
from scipy.optimize import fminbound


## Microphysics constants
g = 9.81 # Gravitational constant, m/s**2
Cp = 1004.0 # Specific heat of dry air at constant pressure, J/kg 
L = 2.5e6 # Latent heat of condensation, J/kg
rho_w = 1e3 # Density of water, kg/m**3
Rd = 287.0 # Gas constant for dry air, J/(kg K)
R = 8.314 # Universal gas constant, J/(mol K)
Mw = 18.0153/1e3 # Molecular weight of water, kg/mol
sigma_w = lambda T: 0.0761 - (1.55e-4)*(T-273.15) # surface tension of water, J/m^2 given T in Kelvin

## NOT CORRECTING FOR NON-CONTINUUM EFFECTS
Dv = 0.3/1e4 # Diffusivity of water vapor in air, m^2/s
#ka = lambda T: 419.*(5.69 + 0.017*(T-273.15))*1e-5 # thermal conductivty of air, W/(m K) given T in Kelvin

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

def Seq(r, r_dry, T, kappa, neg=False):
    '''Equilibrium supersaturation predicted by Kohler theory
    
    Includes optional switch `neg` for inverting the function - useful for
    finding the maxima numerically
    
    following Petters and Kredenweis, 2007
    '''
    A = (2.*Mw*sigma_w(T))/(R*T*rho_w*r)
    B = (r**3 - (r_dry**3))/(r**3 - (r_dry**3)*(1.-kappa))
    if neg:
        return 1.0 - np.exp(A)*B
    else:
        return np.exp(A)*B - 1.0

def kohler_crit(T, r_dry, kappa):
    '''Numerically find the critical radius predicted by kappa Kohler theory'''
    out = fminbound(Seq, r_dry, r_dry*1e3, args=(r_dry, T, kappa, True),
                    xtol=1e-10, full_output=True, disp=0)
    r_crit, s_crit = out[:2]
    s_crit *= -1.0
    return r_crit, s_crit