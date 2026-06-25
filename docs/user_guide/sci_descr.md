# Scientific Description

The simplest tools for describing the growth and evolution of a cloud droplet spectrum
from a given aerosol population are zero-dimensional adiabatic cloud parcel models. By
employing a detailed description of the condensation of ambient water vapor onto growing
droplets, these models accurately describe the activation of a subset of the aerosol
population by predicting how the aerosols modify the maximum supersaturation achieved as
the parcel rises. These models also serve as the theoretical basis for parameterizations
of droplet activation included in modern general circulation models [[Ghan2011](#Ghan2011)].

## Model Formulation

The adiabatic cloud parcel model implemented here is based on models described in the
literature [[Nenes2001](#Nenes2001), [SP2006](#SP2006)] with modifications and improvements; for a full
description see [[Rothenberg2016](#Rothenberg2016)]. The conservation of heat in a parcel rising at constant
velocity $V$ without entrainment is:

$$
\frac{dT}{dt} = -\frac{gV}{c_p} - \frac{L}{c_p}\frac{dw_v}{dt}
$$

where $T$ is the parcel air temperature. Assuming adiabaticity and neglecting entrainment
is appropriate near cloud base, where the majority of droplet activation occurs. Mass
conservation between vapor and liquid phases requires:

$$
\frac{dw_v}{dt} = -\frac{dw_c}{dt}
$$

where $w_v$ and $w_c$ are the mass mixing ratios of water vapor and condensed liquid.
The total condensation rate onto a population of $N_i$ droplets of radius $r_i$
($i = 1, \ldots, n$) is:

$$
\frac{dw_c}{dt} = \frac{4\pi\rho_w}{\rho_a}\sum_{i=1}^n N_i r_i^2 \frac{dr_i}{dt}
$$

The individual droplet growth rate is:

$$
\frac{dr_i}{dt} = \frac{G}{r_i}(S - S_{\text{eq}})
$$

where $G$ is a growth coefficient depending on the physical and chemical properties of
the particle:

$$
G = \left(\frac{\rho_w R T}{e_s D'_v M_w} + \frac{L\rho_w\left[\frac{LM_w}{RT} - 1\right]}{k'_a T}\right)^{-1}
$$

Droplet growth is modulated by the difference between the environmental supersaturation
$S$ and the droplet equilibrium supersaturation $S_{\text{eq}}$ from Köhler theory. The
$\kappa$-Köhler parameterization [[PK2007](#PK2007)] describes hygroscopicity with a single
parameter $\kappa$ relating the water activity of an aqueous solution to the ratio of
solute to water volumes:

$$
\frac{1}{a_w} = 1 + \kappa\frac{V_s}{V_w}
$$

The full $\kappa$-Köhler equation for the equilibrium supersaturation is:

$$
S_{\text{eq}} = \frac{r_i^3 - r_{d,i}^3}{r_i^3 - r_{d,i}^3(1 - \kappa_i)}
\exp\!\left(\frac{2M_w\sigma_w}{RT\rho_w r_i}\right) - 1
$$

where $r_d$ and $r$ are the dry and total (wetted) radii, and the surface tension of
water $\sigma_w = 0.0761 - 1.55 \times 10^{-4}(T - 273.15)$ J m$^{-2}$.

Both the water vapor diffusivity and the thermal conductivity of air are modified for
non-continuum effects as droplets grow:

$$
D'_v = D_v \bigg/ \left(1 + \frac{D_v}{a_c r}\sqrt{\frac{2\pi M_w}{RT}}\right)
$$

$$
k'_a = k_a \bigg/ \left(1 + \frac{k_a}{a_T r \rho_a c_p}\sqrt{\frac{2\pi M_a}{RT}}\right)
$$

where the thermal accommodation coefficient $a_T = 0.96$ and the condensation
coefficient $a_c = 1.0$. Under the adiabatic assumption, the supersaturation evolves
as a balance between condensational heating and adiabatic cooling:

$$
\frac{dS}{dt} = \alpha V - \gamma\frac{dw_c}{dt}
$$

where $\alpha$ and $\gamma$ are weakly temperature- and pressure-dependent:

$$
\alpha = \frac{gM_w L}{c_p R T^2} - \frac{gM_a}{RT}, \qquad
\gamma = \frac{P M_a}{e_s M_w} + \frac{M_w L^2}{c_p R T^2}
$$

The parcel pressure evolves hydrostatically (using virtual temperature to account for
moisture):

$$
\frac{dP}{dt} = \frac{-g P V}{R_d T_v}
$$

These five equations — for $P$, $S$, $w_c$, $w_v$, and $T$ — together with one
radius equation per aerosol bin form a closed system of $7 + n_r$ ODEs that is
integrated forward in time.

## Aerosol Population Specification

The model accepts any population of aerosols representable as a sectional distribution.
The most common representation is a lognormal size distribution:

$$
n_N(r) = \frac{dN}{d\ln r} = \frac{N_t}{\sqrt{2\pi}\ln\sigma_g}
\exp\!\left(-\frac{\ln^2(r/\mu_g)}{2\ln^2\sigma_g}\right)
$$

parameterized by total number concentration $N_t$, geometric mean radius $\mu_g$, and
geometric standard deviation $\sigma_g$. Complex multi-modal distributions are
represented as sums of lognormal modes.

The model discretizes each continuous distribution into $n$ size bins. If no explicit
bounds are specified, $n$ log-spaced bins are placed over $[\mu_g / 10\sigma_g,\,
10\sigma_g \mu_g]$. Typical runs use 100–200 bins per mode; the model shows little
sensitivity to the bin density [[Rothenberg2016](#Rothenberg2016)].

A single $\kappa$ value is prescribed per mode, but the model tracks hygroscopicity
per bin, allowing size-dependent composition and external mixing states.

## Initial Conditions

Before integration begins, the model equilibrates the initial wet-radius spectrum by
numerically solving the Köhler equation (or its linear approximation
$S \approx A/r - \kappa r_d^3/r^3$) for each bin at the initial supersaturation $S_0$.
The resulting equilibrium radii set the initial droplet sizes and determine the initial
liquid water content.

## Numerical Implementation

The ODE system is integrated using the [diffrax](https://github.com/patrick-kidger/diffrax)
library's `Kvaerno5` solver — a fifth-order ESDIRK implicit Runge–Kutta method with
adaptive step control. The right-hand side is compiled via `jax.jit` (XLA) at the first
call; subsequent calls with the same array shapes carry no Python overhead and dispatch
to GPU when requested. A detailed discussion of the solver configuration, tolerances, and
adjoint differentiation is in the [Numerical Methods](numerical_methods.md) guide.

## References

<a id="Nenes2001"></a>**[Nenes2001]** Nenes, A., Ghan, S., Abdul-Razzak, H., Chuang, P. Y., & Seinfeld, J. H.
(2001). Kinetic limitations on cloud droplet formation and impact on cloud albedo.
*Tellus B*, **53**(2), 133–149. doi:[10.1034/j.1600-0889.2001.d01-12.x](https://doi.org/10.1034/j.1600-0889.2001.d01-12.x)

<a id="SP2006"></a>**[SP2006]** Seinfeld, J. H., & Pandis, S. N. (2006). *Atmospheric Chemistry and Physics:
From Air Pollution to Climate Change* (2nd ed.). Wiley.

<a id="Rothenberg2016"></a>**[Rothenberg2016]** Rothenberg, D., & Wang, C. (2016). Metamodeling of droplet activation
for global climate models. *J. Atmos. Sci.*, **73**(4), 1255–1272.
doi:[10.1175/JAS-D-15-0223.1](https://doi.org/10.1175/JAS-D-15-0223.1)

<a id="PK2007"></a>**[PK2007]** Petters, M. D., & Kreidenweis, S. M. (2007). A single parameter
representation of hygroscopic growth and cloud condensation nucleus activity.
*Atmos. Chem. Phys.*, **7**(8), 1961–1971. doi:[10.5194/acp-7-1961-2007](https://doi.org/10.5194/acp-7-1961-2007)

<a id="Ghan2011"></a>**[Ghan2011]** Ghan, S. J., et al. (2011). Droplet nucleation: Physically-based
parameterizations and comparative evaluation. *J. Adv. Model. Earth Syst.*,
**3**(4), M10001. doi:[10.1029/2011MS000074](https://doi.org/10.1029/2011MS000074)
