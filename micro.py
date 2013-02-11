"""
.. module:: parcel
    :synopsis: Useful and commonly used microphysics constants and functions.

.. moduleauthor:: Daniel Rothenberg <darothen@mit.edu>

"""
import numpy as np
from scipy.optimize import fminbound
from scipy.special import erfc, erf


## Microphysics constants
g = 9.81 #: Gravitational constant, m/s**2
Cp = 1004.0 #: Specific heat of dry air at constant pressure, J/kg
L = 2.5e6 #: Latent heat of condensation, J/kg
rho_w = 1e3 #: Density of water, kg/m**3
Rd = 287.0 #: Gas constant for dry air, J/(kg K)
R = 8.314 #: Universal gas constant, J/(mol K)
Mw = 18.0153/1e3 #: Molecular weight of water, kg/mol

# Additional constants from the auxiliary file
Ma = 28.9/1e3 # Molecular weight of dry air, kg/mol
#Dv = 3.e-5 # Diffusivity of water vapor in air, m^2/s
ac = 1.0 # condensation constant
Ka = 2.e-2 # Thermal conductivity of air, J/m/s/K
at = 0.96 # thermal accomodation coefficient
epsilon = 0.622 # molecular weight of water / molecular weight of dry air

## NOT CORRECTING FOR NON-CONTINUUM EFFECTS
Dv = 0.3/1e4 # Diffusivity of water vapor in air, m^2/s
ka = lambda T: 419.*(5.69 + 0.017*(T-273.15))*1e-5 # thermal conductivty of air, W/(m K) given T in Kelvin

## SNIPPETS:
#parcel.Tv = (1. + 0.61*parcel.wv)*parcel['T']
#parcel.rho = parcel.P/(Rd*parcel.Tv)

## AUXILIARY FUNCTIONS
def activation(V, T, P, aerosols):

    ## Originally from Abdul-Razzak 198 w/ Ma. Need kappa formulation
    alpha = (g*Mw*L)/(Cp*R*(T**2)) - (g*Ma)/(R*T)
    gamma = (R*T)/(es(T-273.15)*Mw) + (Mw*(L**2))/(Cp*Ma*T*P)

    ## Condensation effects
    G_a = (rho_w*R*T)/(es(T-273.15)*Dv*Mw)
    G_b = (L*rho_w*((L*Mw/(R*T))-1))/(ka(T)*T)
    G = 1./(G_a + G_b)

    Smis = []
    Sparts = []
    for aerosol in aerosols:
        sigma = aerosol.sigma
        am = aerosol.mu*1e-6
        N = aerosol.N*1e6
        kappa = aerosol.kappa

        fi = 0.5*np.exp(2.5*(np.log(sigma)**2))
        gi = 1.0 + 0.25*np.log(sigma)


        A = (2.*sigma_w(T)*Mw)/(rho_w*R*T)
        _, Smi2 = kohler_crit(T, am, kappa)

        zeta = (2./3.)*A*(np.sqrt(alpha*V/G))
        etai = ((alpha*V/G)**(3./2.))/(N*gamma*rho_w*2.*np.pi)

        ##
        Spa = fi*((zeta/etai)**(1.5))
        Spb = gi*(((Smi2**2)/(etai + 3.*zeta))**(0.75))
        S_part = (1./(Smi2**2))*(Spa + Spb)

        Smis.append(Smi2)
        Sparts.append(S_part)

    Smax = 1./np.sqrt(np.sum(Sparts))

    act_fracs = []
    for Smi, aerosol in zip(Smis, aerosols):
        ui = 2.*np.log(Smi/Smax)/(3.*np.sqrt(2.)*np.log(aerosol.sigma))
        N_act = 0.5*aerosol.N*erfc(ui)
        act_fracs.append(N_act/aerosol.N)

    return Smax, act_fracs

'''
    print Smax
    print "osmotic", N_act, N_act/N
    print "kappa", N_act2, N_act2/N

    x = np.linspace(-3, 3, 1000)
    plot(x, erfc(x))
    plot(u, erfc(u), marker='o', markersize=14)
    plot(u2, erfc(u2), marker='d', color='r', markersize=14)
'''


def sigma_w(T):
    """Calculate the surface tension of water for a given temperature.

    .. math::
        \sigma_w = 0.0761 - 1.55\\times 10^{-4}(T - 273.15)

    :param T: temperature (deg C)
    :type T: float
    :returns: :math:`\sigma_w(T)` (J/m**2)

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

    """
    A = (2.*Mw*sigma_w(T))/(R*T*rho_w*r)
    if kappa == 0.0:
        B = 1.
    else:
        B = (r**3 - (r_dry**3))/(r**3 - (r_dry**3)*(1.-kappa))
    #print A, B, np.exp(A), np.exp(A)*B
    if neg:
        return 1.0 - np.exp(A)*B
    else:
        return np.exp(A)*B - 1.0

def critical_curve(T, r_a, r_b, kappa):
    """Calculates curves of critical wet radius and supersaturations for aerosol
    particles.

    Calls :func:`kohler_crit` for values of ``r_dry`` between ``r_a`` and ``r_b``
    to calculate how the critical supersaturation changes with the dry radius for a
    particle of specified :math:`\\kappa = ` ``kappa``.

    :param T: environment temperature (K)
    :type T: float

    :param r_a, r_b: parcel dry radii left and right bounds (m)
    :type r_a, r_b: float, float

    :param kappa: particle hygroscopicity
    :type kappa: float > 0.0

    :returns: arrays of dry radii, critical radii, and supersaturations

    """
    rs = np.logspace(np.log10(r_a), np.log10(r_b), 200)
    ss = np.array(map(lambda r_dry: kohler_crit(T, r_dry, kappa), rs))
    return rs, ss[:, 0], ss[:, 1]

def kohler_crit(T, r_dry, kappa):
    """Calculates the critical radius and supersaturation of an aerosol particle.

    Passes the negative or inverted supersaturation function from :func:`Seq` to
    a minimum-value optimization routine from SciPy to efficiently calculate the
    radius at which the curve achieves its maximum supersaturation, and that
    supersaturation value.

    :param T: environment temperature (K)
    :type T: float

    :param r_dry: parcel dry radius (m)
    :type r_dry: float

    :param kappa: particle hygroscopicity
    :type kappa: float

    :returns: r_crit - critical radius (value of `r` maximizing :func:`Seq`)
    :returns: s_crit - :func:`Seq` evaluated at ``r_crit``

    """
    '''Numerically find the critical radius predicted by kappa Kohler theory'''
    out = fminbound(Seq, r_dry, r_dry*1e3, args=(r_dry, T, kappa, True),
                    xtol=1e-10, full_output=True, disp=0)
    r_crit, s_crit = out[:2]
    s_crit *= -1.0
    return r_crit, s_crit

def act_fraction(Smax, T, rs, kappa, r_drys, Nis):
    """Calculates the equilibrium activated fraction given the details of a population
    of aerosol sizes.

    NOTE - This works for a *single mode*. In order to study the entire aerosol
    population activated in the parcel model, this will need to be called however
    many times there are modes for each separate mode.
    """
    r_crits, s_crits = zip(*[kohler_crit(T, r_dry, kappa) for r_dry in r_drys])
    s_crits = np.array(s_crits)
    r_crits = np.array(r_crits)

    activated_eq = (Smax >= s_crits)
    activated_kn = (rs >= r_crits)

    N_tot = np.sum(Nis)

    eq_frac = np.sum(Nis[activated_eq])/N_tot
    kn_frac = np.sum(Nis[activated_kn])/N_tot

    return eq_frac, kn_frac

def r_eff(rho, wc, Ni):
    """Calculates the cloud droplet effective radius given the parcel liquid
    water mixing ratio and the number of activated droplets, as well as the parcel
    air density

    Assuming the droplet population is monodisperse or close to it, the cloud droplet
    effective radius can be computed by

    .. math::
        r_{\\text{eff}} = \left(\\frac{3 \\rho_a w_c}{4 \pi N_i \\rho_w}\\right)^{1/3}

    :param rho: parcel air density (kg/m^3)
    :type rho: float

    :param wc: liquid water mixing ratio (kg/kg)
    :type wc: float

    :param Ni: droplet number concentration (1/m^3)
    :type Ni: float

    :returns: droplet effective radius (m^3)
    """
    return (3.*rho*wc/(4.*np.pi*rho_w*Ni))**(1./3.)