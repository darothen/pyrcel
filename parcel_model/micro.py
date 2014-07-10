"""
.. module:: parcel
    :synopsis: Useful and commonly used microphysics/thermodynamics functions.

.. moduleauthor:: Daniel Rothenberg <darothen@mit.edu>

"""
import numpy as np
from scipy.optimize import fminbound, bisect
from scipy.special import erfc, erf

from constants import *

##########################
## THERMODYNAMIC FUNCTIONS

def dv_cont(T, P):
    """Diffusivity of water vapor in air, neglectiving non-continuum effects.

    See :func:`dv` for equation.

    **Args**:
        * *T* -- ambient temperature of air surrounding droplets, K
        * *P* -- ambient pressure of surrouding air, Pa

    **Returns**:
        :math:`D_v(T, P)` in m^2/s

    """
    P_atm = P*1.01325e-5 # Pa -> atm
    return 1e-4*(0.211/P_atm)*((T/273.)**1.94)

def dv(T, r, P):
    """Diffusivity of water vapor in air, modified for non-continuum effects.

    The diffusivity of water vapor in air as a function of temperature and pressure
    is given by

    .. math::
        \\begin{equation}
        D_v = 10^{-4}\\frac{0.211}{P}\left(\\frac{T}{273}\\right)^{1.94} \\tag{SP2006, 17.61}
        \end{equation}

    where :math:`P` is in atm. Accounting for continuum regime effects which can
    influence the evolving aerosol and droplets produces a modified diffusivity
    implemented here

    .. math::
        \\begin{equation}
        D'_v = \\frac{D_v}{1+ \\frac{D_v}{\\alpha_c r_p}\left(\\frac{2\pi M_w}{RT}\\right)^{1/2}} \\tag{SP2006, 17.62}
        \end{equation}

    where :math:`\\alpha_c` = ``ac`` and :math:`r_p` = ``r``.

    **Args**:
        * *T* -- ambient temperature of air surrounding droplets, K
        * *r* -- radius of aerosol/droplet, m
        * *P* -- ambient pressure of surrounding air, Pa

    **Returns**:
        :math:`D'_v(T, r, P)` in m^2/s

    """
    #P_atm = P*1.01325e-5 # Pa -> atm
    #dv_t = 1e-4*(0.211/P_atm)*((T/273.)**1.94)
    dv_t = dv_cont(T, P)
    denom = 1.0 + (dv_t/(ac*r))*np.sqrt((2.*np.pi*Mw)/(R*T))
    return dv_t/denom

def rho_air(T, P, RH=1.0):
    """ Calculate the density of moist air by assuming a given relative humidty,
    temperature, and pressure.

    Uses the traditional formula (3.41) from Petty.

    .. math::
        \\begin{equation}
        \rho_{\text{air}} = \frac{P}{R_dT_v}
        \end{equation}

    **Args**:
        * *T* -- ambient air temperature, K
        * *P* -- ambient air pressure, Pa
        * *RH* -- relative humidity, decimal

    **Returns**:
        density of air in kg/m^3

    """
    qsat = RH*0.622*(es(T-273.15)/P)
    Tv = T*(1.0 + 0.61*qsat)
    rho_a = P/Rd/Tv # air density
    
    return rho_a

def es(T_c):
    """Calculates the saturation vapor pressure over water for a given temperature.

    Uses the Bolton (1980) empirical fit, which is accurate to :math:`0.1\%` over the
    temperature range :math:`-30^oC \leq T \leq 35^oC`,

    .. math::
        \\begin{equation}
        e_s(T) = 611.2 \exp\left(\\frac{17.67T}{T + 243.5}\\right) \\tag{RY1989, 2.17}
        \end{equation}

    where :math:`e_s` is in Pa and :math:`T` is in degrees C.

    **Args**:
        * *T_c* -- ambient air temperature, degrees C

    **Returns**:
        :math:`e_s(T)` in Pa

    """
    return 611.2*np.exp(17.67*T_c/(T_c+243.5))

def ka_cont(T):
    """Thermal conductivity of air, neglecting non-continuum effects.

    See :func:`ka` for equation.

    **Args**:
        * *T* -- ambient air temperature surrounding droplet, K

    **Returns**:
        :math:`k_a(T)` in J/m/s/K

    """
    return 1e-3*(4.39 + 0.071*T)

def ka(T, rho, r):
    """Thermal conductivity of air, modified for non-continuum effects.

    The thermal conductivity of air is given by

    .. math::
        \\begin{equation}
        k_a = 10^{-3}(4.39 + 0.071T) \\tag{SP2006, 17.71}
        \end{equation}

    Modified to account for non-continuum effects yields the equation

    .. math::
        \\begin{equation}
        k'_a = \\frac{k_a}{1 + \\frac{k_a}{\\alpha_t r_p \\rho C_p} \\frac{2\pi M_a}{RT}^{1/2}} \\tag{SP2006, 17.72}
        \end{equation}

    where :math:`\\alpha_t` is a thermal accommodation coefficient.

    **Args**:
        * *T* -- ambient air temperature, K
        * *rho* -- ambient air density, kg/m^3
        * *r* -- particle/droplet radius, m

    **Returns**:
        :math:`k'_a(T, \\rho, r_p)` in J/m/s/K

    """
    ka_t = ka_cont(T)
    denom = 1.0 + (ka_t/(at*r*rho*Cp))*np.sqrt((2.*np.pi*Ma)/(R*T))
    return ka_t/denom

def sigma_w(T):
    """Calculates the surface tension of water for a given temperature.

    .. math::
        \\begin{equation}
        \sigma_w = 0.0761 - 1.55\\times 10^{-4}(T - 273.15)
        \end{equation}

    **Args**:
        * *T* -- ambient air temperature, degrees K

    **Returns**:
        :math:`\sigma_w(T)` in J/m^2

    """
    return 0.0761 - (1.55e-4)*(T-273.15)

##########################
## KOHLER THEORY FUNCTIONS

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
    this method will return the supersaturation as a decimal with respect to 1.0,

    .. math::
        S_\\text{eq} = S - 1.0

    **Args**:
        * *r* -- droplet radius, m
        * *r_dry* -- dry particle radius, m
        * *T* -- ambient air temperature, K
        * *kappa* -- particle hygroscopicity parameter
        * *neg* -- (optional) boolean flag indicating whether to return the negative of the calculation

    **Returns**:
        :math:`S_\\text{eq}`

    """
    A = (2.*Mw*sigma_w(T))/(R*T*rho_w*r)
    if kappa == 0.0:
        B = 1.
    else:
        B = (r**3 - (r_dry**3))/(r**3 - (r_dry**3)*(1.-kappa))
    if neg:
        return 1.0 - np.exp(A)*B
    else:
        return np.exp(A)*B - 1.0

def kohler_crit(T, r_dry, kappa):
    """Calculates the critical radius and supersaturation of an aerosol particle.

    Passes the negative or inverted supersaturation function from :func:`Seq` to
    a minimum-value optimization routine from SciPy to efficiently calculate the
    radius at which the curve achieves its maximum supersaturation, and that
    supersaturation value.

    **Args**:
        * *T* -- ambient air temperature, K
        * *r_dry* -- dry particle radius, m
        * *kappa* -- particle hygroscopicity parameter

    **Returns**:
        Tuple of :math:`(r_\\text{crit},\, S_\\text{crit})`, the critical radius (m)
        and supersaturation of the aerosol droplet.

    """
    out = fminbound(Seq, r_dry, r_dry*1e4, args=(r_dry, T, kappa, True),
                    xtol=1e-10, full_output=True, disp=0)
    r_crit, s_crit = out[:2]
    s_crit *= -1.0 # multiply by -1 to undo negative flag for Seq
    return r_crit, s_crit


def kappa_kohler_crit(T, r_dry, kappa):
    """Calculates the critical radius and supersaturation of an aerosol particle.

    Similar to :func:`kohler_crit`, but avoids a numerical iteration by using the
    approximation

    .. math::
        \\begin{equation}
        S_\\text{eq} = \\frac{A}{r} - \kappa\\frac{r^3_\\text{dry}{r^3}
        \end{equation}

    where :math:`A = 2M_w\sigma_w(T)/(R T \\rho_w). The critical values are set by
    taking :math:`dS_\\text{eq}/dr = 0` and solving for :math:`r`.

    **Args**:
        * *T* -- ambient air temperature, K
        * *r_dry* -- dry particle radius, m
        * *kappa* -- particle hygroscopicity parameter

    **Returns**:
        Tuple of :math:`(r_\\text{crit},\, S_\\text{crit})`, the critical radius (m)
        and supersaturation of the aerosol droplet.

    .. note::

        It's almost always preferable to use :func:`kohler_crit` as it will yield the
        most accurate calculation for the critical supersaturation and radius.

    """
    A = 4.*Mw*sigma_w(T)/R/T/rho_w
    rd3 = r_dry**3
    r_crit = np.sqrt(3.*kappa*rd3/A)
    s_crit = np.sqrt(4.*(A**3)/(27.*kappa*rd3))
    return r_crit, s_crit

def critical_curve(T, r_a, r_b, kappa):
    """Calculates curves of critical wet radius and supersaturations for aerosol
    particles.

    Calls :func:`kohler_crit` for values of ``r_dry`` between ``r_a`` and ``r_b``
    to calculate how the critical supersaturation changes with the dry radius for a
    particle of specified :math:`\\kappa = ` ``kappa``.

    **Args**:
        * *T* -- ambient air temperature, K
        * *r_a*, *r_b* -- left/right bounds of parcel dry radii, m
        * *kappa* -- particle hygroscopicity parameter

    **Returns**:
        Three arrays,
            * particle dry radii between ``r_a`` and ``r_b``,
            * critical radii for those particles
            * critical supersaturation

    """
    rs = np.logspace(np.log10(r_a), np.log10(r_b), 200)
    ss = np.array(map(lambda r_dry: kohler_crit(T, r_dry, kappa), rs))
    return rs, ss[:, 0], ss[:, 1]

###############
## MICROPHYSICS

def r_eff(rho, wc, Ni):
    """Calculates the cloud droplet effective radius given the parcel liquid
    water mixing ratio and the number of activated droplets, as well as the parcel
    air density

    Assuming the droplet population is monodisperse or close to it, the cloud droplet
    effective radius can be computed by

    .. math::
        \\begin{equation}
        r_{\\text{eff}} = \left(\\frac{3 \\rho_a w_c}{4 \pi N_i \\rho_w}\\right)^{1/3}
        \end{equation}

    **Args**:
        * *rho* -- parcel air density, kg/m^3
        * *wc* -- liquid water mixing ratio, kg/kg
        * *Ni* -- droplet number concentration, m^-3

    **Returns**:
        Cloud droplet effective radius, m

    .. warning::

        Not completely implemented yet.

    """
    return (3.*rho*wc/(4.*np.pi*rho_w*Ni))**(1./3.)
