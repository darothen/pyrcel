.. _sci_descr:

Scientific Description
======================

The simplest tools available for describing the growth and evolution of
a cloud droplet spectrum from a given population of aerosols are based
on zero-dimensional, adiabatic cloud parcel models. By employing a
detailed description of the condensation of ambient water vapor onto the
growing droplets, these models can accurately describe the activation of
a subset of the aerosol population by predicting how the presence of the
aerosols in the updraft modify the maximum supersaturation achieved as
the parcel rises. Furthermore, these models serve as the theoretical
basis or reference for parameterizations of droplet activation which are
included in modern general circulation models ([Ghan2011]_) .

The complexity of these models varies with the range of physical processes
one wishes to study. At the most complex end of the spectrum, one might wish
to accurately resolve chemical transfer between the gas and aqueous phase in
addition to physical transformations such as collision/coalescence. One could
also add ice-phase processes to such a model.

Model Formulation
-----------------

The adiabatic cloud parcel model implemented here is based on
models described in the literature ([Nenes2001]_, [SP2006]_,) with some modifications and improvements. For a full description of the parcel model, please see ([Rothenberg2016]_)
The conservation of heat in a parcel of air rising at constant
velocity :math:`V` without entrainment can be written as

.. math:: \frac{dT}{dt} = -\frac{gV}{c_p} - \frac{L}{c_p}\frac{d w_v}{dt}
    :label: dTdt

where :math:`T` is the parcel air temperature. Assuming adiabaticity
and neglecting entrainment is suitable for studying cloud droplet
formation near the cloud base, where the majority of droplet activation
occurs. Because the mass of water must be conserved as it changes from
the vapor to liquid phase, the relationship

.. math:: \frac{d w_v}{dt} = - \frac{dw_c}{dt}
    :label: dw_vdt

must hold, where :math:`w_v` and :math:`w_c` are the mass mixing ratios
of water vapor and liquid water (condensed in droplets) in the parcel.
The rate of change of water in the liquid phase in the parcel is
governed solely by condensation onto the existing droplet population.
For a population of :math:`N_i` droplets of radius :math:`r_i`, where
:math:`i=1,\dots,n`, the total condensation rate is given by

.. math:: \frac{dw_c}{dt} = \frac{4\pi \rho_w}{\rho_a}\sum\limits_{i=1}^nN_ir_i^2\frac{dr_i}{dt}
    :label: dw_cdt

Here, the particle growth rate, :math:`\frac{dr_i}{dt}` is calculated as

.. math:: \frac{dr_i}{dt} = \frac{G}{r_i}(S-S_{eq})
    :label: dr_idt

where :math:`G` is a growth coefficient which is a function of the
physical and chemical properties of the particle receiving condensate,
given by

.. math:: G = \left(\frac{\rho_w R T}{e_s D'_v M_w} + \frac{L\rho_w[(LM_w/RT) - 1]}{k'_a T}\right)^{-1}
    :label: G

Droplet growth via condensation is modulated by the difference between
the environmental supersaturation, :math:`S`, and the droplet
equilibrium supersaturation, :math:`S_{eq}`, predicted from Kohler
theory. To account for differences in aerosol chemical properties which
could affect the ability for particles to uptake water, the
:math:`\kappa`-Köhler theory parameterization ([PK2007]_) is employed in the
model. :math:`\kappa`-Kohler theory utilizes a single parameter to
describe aerosol hygroscopicity, and is widely employed in modeling of
aerosol processes. The hygroscopicity parameter :math:`\kappa` is
related to the water activity of an aqueous aerosol solution by

.. math:: \frac{1}{a_w} = 1 + \kappa\frac{V_s}{V_w}

where :math:`V_s` and :math:`V_w` are the volumes of dy particulate
matter and water in the aerosol solution. With this parameter, the full
:math:`\kappa`-Kohler theory may be expressed as

.. math:: S_{eq} = \frac{r_i^3 - r_{d,i}^3}{r_i^3 - r_{d,i}^3(1-\kappa_i)}\exp\left( \frac{2M_w\sigma_w}{RT\rho_w r_i} \right) - 1
    :label: S_eq

where :math:`r_d` and :math:`r` are the dry aerosol particle size and
the total radius of the wetted aerosol. The surface tension of water,
:math:`\sigma_w`, is dependent on the temperature of the parcel such
that :math:`\sigma_w = 0.0761 - 1.55\times 10^{-4}(T-273.15)`
J/m\ :math:`^2` . Both the diffusivity and thermal conductivity of air
have been modified in the growth coefficient equation to account for
non-continuum effects as droplets grow, and are given by the expressions

.. math:: D'_v = D_v\bigg/\left(1 + \frac{D_v}{a_c r}\sqrt{\frac{2\pi M_w}{RT}}\right)

and

.. math:: k'_a = k_a\bigg/\left(1 + \frac{k_a}{a_T r \rho_a c_p}\sqrt{\frac{2\pi M_a}{RT}} \right)

In these expressions, the thermal accommodation coefficient,
:math:`a_T`, is assumed to be :math:`0.96` and the condensation
coefficient, :math:`a_c` is taken as unity (see :ref:`Constants <constants>`).
Under the adiabatic assumption, the evolution of the parcel’s
supersaturation is governed by the balance between condensational
heating as water vapor condenses onto droplets and cooling induced by
the parcel’s vertical motion,

.. math:: \frac{dS}{dt} = \alpha V - \gamma\frac{w_c}{dt}
    :label: dSdt

where :math:`\alpha` and :math:`\gamma` are functions which are weakly
dependent on temperature and pressure :

.. math:: \alpha = \frac{gM_wL}{c_pRT^2} - \frac{gM_a}{RT}

.. math:: \gamma = \frac{PM_a}{e_sM_w} + \frac{M_wL^2}{c_pRT^2}

The parcel’s pressure is predicted using the hydrostatic relationship,
accounting for moisture by using virtual temperature (which can always
be diagnosed as the model tracks the specific humidity through the mass
mixing ratio of water vapor),

.. math:: \frac{dP}{dt} = \frac{-g P V}{R_d T_v}
    :label: dPdt

The equations :eq:`dPdt`, :eq:`dSdt`, :eq:`dw_cdt`, :eq:`dw_vdt`,
and :eq:`dTdt` provide a simple, closed system of ordinary
differential equations which can be numerically integrated forward in
time. Furthermore, this model formulation allows the simulation of an
arbitrary configuration of initial aerosols, in terms of size, number
concentration, and hygroscopicity. Adding additional aerosol size bins
is simply accomplished by tracking one additional size bin in the system
of ODE’s. The important application of this feature is that the model
can be configured to simulate both internal or external mixtures of
aerosols, or some combination thereof.

Model Implementation and Procedure
----------------------------------

The parcel model described in the previous section was implemented using
a modern modular and object-oriented software engineering framework.
This framework allows the model to be simply configured with myriad
initial conditions and aerosol populations. It also enables model
components - such as the numerical solver or condensation
parameterization - to be swapped and replaced. Most importantly, the use
of object-oriented techniques allows the model to be incorporated into
frameworks which grossly accelerate the speed at which the model can be
evaluated. For instance, although models like the one developed here are
relatively cheap to execute, large ensembles of model runs have been
limited in scope to several hundred or a thousand runs. However, the
framework of this particular parcel model implementation was designed
such that it could be run as a black box as part of a massively-parallel
ensemble driver.

To run the model, a set of initial conditions needs to be specified,
which includes the updraft speed, the parcel’s initial temperature,
pressure, and supersaturation, and the aerosol population. Given these
parameters, the model calculates an initial equilibrium droplet spectrum
by computing the equilibrium wet radii of each aerosol. This is calculated
numerically from the Kohler equation for each aerosol/proto-droplet, or
numerically by employing the typical Kohler theory approximation

.. math:: S \approx \frac{A}{r} - \kappa\frac{r_d^3}{r^3}

These wet radii are used as the initial droplet radii in the simulation.

Once the initial conditions have been configured, the model is
integrated forward in time with a numerical solver (see :func:`ParcelModel.run`
for more details). The available solvers wrapped here are:

- LSODA(R)
- LSODE
- (C)VODE

During the model integration, the size representing each aerosol bin is
allowed to grow via condensation, producing something akin to a moving
grid.  In the future, a fixed Eulerian
grid will likely be implemented in the model for comparison.


Aerosol Population Specification
--------------------------------

The model may be supplied with any arbitrary population of aerosols,
providing the population can be approximated with a sectional
representation. Most commonly, aerosol size distributions are
represented with a continuous lognormal distribution,

.. math:: n_N(r) = \frac{dN}{d \ln r} = \frac{N_t}{\sqrt{2\pi}\ln \sigma_g}\exp\left(-\frac{ \ln^2(r/\mu_g)}{2\ln^2\sigma_g}\right)
    :label: lognormal

which can be summarized with the set of three parameters,
:math:`(N_t, \mu_g, \sigma_g)` and correspond, respectively, to the
total aerosol number concentration, the geometric mean or number mode
radius, and the geometric standard deviation. Complicated multi-modal
aerosol distributions can often be represented as the sum of several
lognormal distributions. Since the parcel model describes the evolution
of a discrete aerosol size spectrum, can be broken into :math:`n` bins,
and the continuous aerosol size distribution approximated by taking the
number concentration and size at the geometric mean value in each bin,
such that the discrete approximation to the aerosol size distribution
becomes

.. math:: n_{N,i}(r_i) = \sum\limits_{i=1}^n\frac{N_i}{\sqrt{2\pi}\ln\sigma_g}\exp\left(-\frac{\ln^2(r_i/\mu_g)}{2\ln^2\sigma_g}\right)

If no bounds on the size range of :math:`r_i` is specified, then the
model pre-computes :math:`n` equally-spaced bins over the logarithm of
:math:`r`, and covers the size range :math:`\mu_g/10\sigma_g` to
:math:`10\sigma_g\mu_g`. It is typical to run the model with :math:`200`
size bins per aerosol mode. Neither this model nor similar ones exhibit
much sensitivity towards the density of the sectional discretization .

Typically, a single value for hygroscopicity, :math:`\kappa` is
prescribed for each aerosol mode. However, the model tracks a
hygroscopicity parameter for each individual size bin, so size-dependent
aerosol composition can be incorporated into the aerosol population.
This representation of the aerosol population is similar to the external
mixing state assumption. An advantage to using this representation is
that complex mixing states can be represented by adding various size
bins, each with their own number concentration and hygroscopicity.

.. topic:: Reference

References
----------

.. [Nenes2001] Nenes, A., Ghan, S., Abdul-Razzak, H., Chuang, P. Y. & Seinfeld, J. H. Kinetic limitations on cloud droplet formation and impact on cloud albedo. Tellus 53, 133–149 (2001).

.. [SP2006] Seinfeld, J. H. & Pandis, S. N. Atmospheric Chemistry and Physics: From Air Pollution to Climate Change. Atmos. Chem. Phys. 2nd, 1203 (Wiley, 2006).

.. [Rothenberg2016] Daniel Rothenberg and Chien Wang, 2016: Metamodeling of Droplet Activation for Global Climate Models. *J. Atmos. Sci.*, **73**, 1255–1272. doi: http://dx.doi.org/10.1175/JAS-D-15-0223.1

.. [PK2007] Petters, M. D. & Kreidenweis, S. M. A single parameter representation of hygroscopic growth and cloud condensation nucleus activity. Atmos. Chem. Phys. 7, 1961–1971 (2007).

.. [Ghan2011] Ghan, S. J. et al. Droplet nucleation: Physically-based parameterizations and comparative evaluation. J. Adv. Model. Earth Syst. 3, M10001 (2011).