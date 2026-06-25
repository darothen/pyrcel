# pyrcel: cloud parcel model

**pyrcel** is a zero-dimensional adiabatic cloud parcel model for studying aerosol–cloud
interactions. It simulates the condensational growth of an arbitrary aerosol population as
a parcel rises adiabatically, predicting the peak supersaturation and the number of
droplets activated at cloud base.

[![DOI](https://zenodo.org/badge/12927551.svg)](https://zenodo.org/badge/latestdoi/12927551)
[![PyPI](https://badge.fury.io/py/pyrcel.svg)](https://badge.fury.io/py/pyrcel)
[![CI](https://github.com/darothen/pyrcel/actions/workflows/ci.yml/badge.svg)](https://github.com/darothen/pyrcel/actions/workflows/ci.yml)
[![Documentation Status](https://readthedocs.org/projects/pyrcel/badge/?version=stable)](https://pyrcel.readthedocs.io/en/stable/)

![Parcel model simulation: supersaturation and aerosol size evolution](assets/figures/basic_run.png)

!!! warning "Version 2.0 Notice"

    This is **pyrcel v2.0**, a major new release with more features, greater flexibility, and a JAX-based differentiable kernel. However, there are several **breaking changes** compared to version 1.3.x.

    Please review the [migration guide](user_guide/migration.md) to update your code.

    If you wish to continue using the legacy (v1.3.x) model, you can install it from PyPI by pinning the version:

    ```shell
    pip install "pyrcel<2"
    ```


## Highlights

**Differentiable.** Every output is a smooth function of its inputs. Compute
$\partial S_\text{max}/\partial V$, $\partial N_d/\partial \kappa$, or full Jacobians
through the ODE using `jax.grad` and `jax.jacfwd`.

**GPU-capable.** Pass `device="gpu"` to `ParcelModel` and the integration dispatches
to CUDA with no code changes. Particularly useful for ensemble or batch simulations!

**Batchable.** Use `jax.vmap` to map a single compiled kernel over an ensemble of
initial conditions — thousands of parcel trajectories in a single call.

**Pure Python.** No compiled extensions or conda dependencies. Install with `pip`.

## Quick install

```bash
uv add pyrcel          # recommended
pip install pyrcel     # pip alternative
```

See the [Installation guide](getting_started/installation.md) for GPU, editable,
and GitHub installs.

---

## Minimal example

```python
import pyrcel as pm

sulfate = pm.AerosolSpecies(
    "sulfate", pm.Lognorm(mu=0.05, sigma=2.0, N=1000.0), kappa=0.54, bins=50
)
model = pm.ParcelModel([sulfate], V=1.0, T0=283.0, S0=-0.02, P0=85000.0, console=True)
output = model.run(t_end=300.0, output_dt=10.0, terminate=True, live=True)

print(f"S_max  = {output.summary['S_max']*100:.3f} %")
print(f"N_act  = {output.Nd:.3e} m⁻³")
```

??? note "Click to expand sample output"

    ```
    Parcel model (JAX / diffrax)
    ------------------------------------------------------------------------
      Backend   JAX 0.10.1 | cpu:0 | float64=on

      Configuration
        V            1 m/s              accom        1
        T0           283.00 K           P0           850.0 hPa
        S0           -0.02000 (-2.00 %)    bins         50

      species           κ     N [cm⁻³]   bins  size
      --------------------------------------------------------------------
      sulfate       0.540       1002.4     50  μ=0.05 μm  σ=2

      Equilibrated initial state
        equilibration   520.3 ms  |Seq-S0|_max = 1.49e-16
           P    850.0 hPa     T  283.00 K    wv     8.84 g/kg    wc 1.05e-04 g/kg     S -0.02000
    ------------------------------------------------------------------------

      Integration plan
        t_end cap     300 s
        output_dt     10 s
        terminate     yes (+10 m past S_max)
        solver        Kvaerno5 + PIDController  rtol=1e-07  max_steps=100000
        progress      live chunk loop (10 s chunks)
    ------------------------------------------------------------------------

      Integration loop

        step     time  walltime  Δwalltime |     z       T       S
       ------------------------------------|----------------------
           1    0.00s     0.00s      0.00s |   0.0  283.00  -2.00%
           2   10.00s     2.45s      2.45s |  10.0  282.90  -1.53%
           3   20.00s     2.47s      0.02s |  20.0  282.80  -1.05%
           4   30.00s     2.49s      0.02s |  30.0  282.71  -0.58%
           5   40.00s     2.52s      0.03s |  40.0  282.61  -0.12%
           6   50.00s     2.56s      0.04s |  50.0  282.52   0.25%
           7   60.00s     2.58s      0.02s |  60.0  282.47   0.20%
       ---- end of integration loop ----

    [pyrcel] Termination: S_max = 0.2691 % at t = 52.50 s (z = 52.5 m); stopped at t = 62.50 s (z = 62.5 m).

      Simulation summary
    ------------------------------------------------------------------------
      S_max = 0.2691 %  at t = 52.50 s (T = 282.52 K, z = 50.0 m)
           species   eq_act   kn_act      N_act          N
           sulfate    0.637    0.570 638373012.4 1002378640.4
      total activated fraction = 0.637
    ------------------------------------------------------------------------

    S_max  = 0.269 %
    N_act  = 6.384e+08 m⁻³
    ```

---

## Cite

If you use pyrcel for research, please cite the original manuscrip where
the model was detailed:

> Rothenberg, D., & Wang, C. (2016). Metamodeling of droplet activation for global
> climate models. *J. Atmos. Sci.*, **73**(4), 1255–1272.
> doi:[10.1175/JAS-D-15-0223.1](https://doi.org/10.1175/JAS-D-15-0223.1)

Additionally, please consider citing the bespoke DOI for the
[specific release version of pyrcel](https://zenodo.org/records/20693507)
that you used during your research (or the base version you modified). This
allows us to track adoption and use of specific model versions over time.

---

## Contents

<div class="grid cards" markdown>

- :material-rocket-launch: **[Getting Started](getting_started/installation.md)**

    Install pyrcel and run your first simulation in minutes.

- :material-book-open-variant: **[User Guide](user_guide/sci_descr.md)**

    Scientific background, migration from v1, numerical methods, and GPU setup.

- :material-code-braces: **[Examples](examples/basic_run.md)**

    Worked examples with rendered output: basic runs, live integration,
    ensemble sweeps, and gradient-based sensitivity analysis.

- :material-api: **[API Reference](api/model.md)**

    Complete reference for all public classes and functions.

</div>
