"""
.. module:: parcel
    :synopsis: Useful and commonly used microphysics constants and functions.

.. moduleauthor:: Daniel Rothenberg <darothen@mit.edu>

"""
import numpy as np
from scipy.optimize import fminbound

## Microphysics constants
g = 9.81 #: Gravitational constant, m/s**2
Cp = 1004.0 #: Specific heat of dry air at constant pressure, J/kg
L = 2.5e6 #: Latent heat of condensation, J/kg
rho_w = 1e3 #: Density of water, kg/m**3
Rd = 287.0 #: Gas constant for dry air, J/(kg K)
R = 8.314 #: Universal gas constant, J/(mol K)
Mw = 18.0153/1e3 #: Molecular weight of water, kg/mol

## NOT CORRECTING FOR NON-CONTINUUM EFFECTS
#Dv = 0.3/1e4 # Diffusivity of water vapor in air, m^2/s
#ka = lambda T: 419.*(5.69 + 0.017*(T-273.15))*1e-5 # thermal conductivty of air, W/(m K) given T in Kelvin

## AUXILIARY FUNCTIONS
def sigma_w(T):
    """Calculate the surface tension of water for a given temperature.

    .. math::
        \sigma_w = 0.0761 - 1.55\\times 10^{-4}(T - 273.15)

    :param T: temperature (deg C)
    :type T: float
    :returns: :math:`\sigma_w(T)` (J/m**2)
    :rtype: float:

    """
    return 0.0761 - (1.55e-4)*(T-273.15)

def es(T):
    """Calculates the saturation vapor pressure over water for a given temperature.

    Implements Formula 2.17 from Rogers and Yau,

    .. math::

        e_s(T) = 6.112\exp\\frac{17.67T}{T+243.5}

    where *T* is the temperature in degrees Celsius.

    :param T: temperature (deg C)
    :type T: float
    :returns: :math:`e_s(T)` (Pa)
    :rtype: float

    """
    return 611.2*np.exp(17.67*T/(T+243.5))

def Seq(r, r_dry, T, kappa, neg=False):
    """Calculates the equilibrium supersaturation (relative to 100% RH) over an
    aerosol parcel of given dry/wet radius and of specified hygroscopicity in a
    parcel of a specific temperature.

    Following the parameterization by Petters and Kredenweis (2007), classical
    Kohler theory can be modified to parameterize the hygroscopicity of an aerosol
    particle using a parameter, :math:`\kappa`. The modified theory predicts that
    the supersaturation with respect to a given aerosol particle is,

    .. math::
        S_\\text{eq} &= a_w \exp \\left( \\frac{2\sigma_{w} M_w}{RT\\rho_w r} \\right) \\\\
        a_w &= \\left(1 + \kappa\\left(\\frac{r_d}{r}^3\\right) \\right)^{-1}

    with the relevant thermodynamic properties of water defined elsewhere in this
    module, :math:`r_d` is the particle dry radius (``r_dry``), :math:`r` is the
    radius of the droplet containing the particle (``r``), :math:`T` is the temperature
    of the environment (``T``), and :math:`\kappa` is the hygroscopicity parameter
    of the particle (``kappa``).

    This method has been extended to supply the *negative* of the supersaturation if
    specified using the argument ``neg``; this is useful when attempting to numerically
    estimate the particle's critical radius, as done in :func:`kohler_crit`. Otherwise,
    this method will return the supersaturation as a decimal with resepct to 1.0,

    .. math::
        S_\\text{eq} = S - 1.0

    :param r: droplet radius (m)
    :type r: float
    :param r_dry: particle dry radius (m)
    :type r_dry: float
    :param T: environment temperature (K)
    :type T: float
    :param kappa: particle hygroscopicity
    :type kappa: float > 0.0
    :param neg: flag indicating to return the negative of the supersaturation, default ``False``
    :type neg: boolean -- optional

    :returns: :math:`S_\\text{eq}`
    :rtype: float

    """
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