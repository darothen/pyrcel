pyrcel: cloud parcel model
==========================

![sample parcel model run](docs/figs/model_example.png)

[![DOI](https://zenodo.org/badge/12927551.svg)](https://zenodo.org/badge/latestdoi/12927551)[![PyPI Version](https://badge.fury.io/py/pyrcel.svg)](https://badge.fury.io/py/pyrcel)[![CircleCI Build Status](https://circleci.com/gh/darothen/pyrcel/tree/master.svg?style=svg)](https://circleci.com/gh/darothen/pyrcel/tree/master)[![Documentation Status](https://readthedocs.org/projects/pyrcel/badge/?version=stable)](http://pyrcel.readthedocs.io/en/latest/index.html)


This is an implementation of a simple, adiabatic cloud parcel model for use in
aerosol-cloud interaction studies. [Rothenberg and Wang (2016)](http://journals.ametsoc.org/doi/full/10.1175/JAS-D-15-0223.1) discuss the model in detail and its improvements
 and changes over [Nenes et al (2001)][nenes2001]:

* Implementation of κ-Köhler theory for condensation physics ([Petters and
Kreidenweis, 2007)][pk2007]
* Extension of model to handle arbitrary sectional representations of aerosol
populations, based on user-controlled empirical or parameterized size distributions
* Improved, modular numerical framework for integrating the model, including bindings
to several different stiff integrators:
    - ~~`lsoda` - [scipy ODEINT wrapper](http://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.odeint.html)~~
    - ~~`vode, lsode*, lsoda*` - ODEPACK via [odespy][hplgit]~~
    - `cvode` - SUNDIALS via [Assimulo](http://www.jmodelica.org/assimulo_home/index.html#)

among other details. It also includes a library of droplet activation routines and scripts/notebooks for evaluating those schemes against equivalent calculations done with the parcel model.

> [!WARNING]
> As of version 1.3, we no longer support any ODE solver backends other than `cvode`.
> All publications using this model have used this backend, so users shouldn't expect
> any inconsistencies with historical results. A future version is planned to add a new
> suite of ODE solvers from the [diffrax][diffrax] toolkit.

Updated code can be found the project [github repository](https://github.com/darothen/pyrcel). If you'd like to use this code or have any questions about it, please [contact the author][author_email]. In particular, if you use this code for research purposes, be sure to carefully read through the model and ensure that you have tweaked/configured it for your purposes (i.e., modifying the accomodation coefficient); other derived quantities).

[Detailed documentation is available](http://pyrcel.readthedocs.org/en/latest/index.html), including a [scientific description](http://pyrcel.readthedocs.org/en/latest/sci_descr.html), [installation details](http://pyrcel.readthedocs.org/en/latest/install.html), and a [basic example](http://pyrcel.readthedocs.org/en/latest/examples/basic_run.html) which produces a figure like the plot at the top of this page.

Quick Start
-----------

As of February, 2025, we provide an ultra simple way to run `pyrcel` without any installation
or setup using [`pixi`](https://pixi.sh/latest/).
`pixi` is an all-in-one package management tool that makes handling complex environment
setup and dependencies extremely easy.

Clone or download this repo, then **cd** into the top-level folder from a terminal.
From there, execute:

``` shell
$ pixi run run_parcel examples/simple.yml
```

This will automatically prepare an environment with all of `pyrcel`'s dependencies installed,
and then run an example model setup.
The first time the model runs, it may take a few second after invoking the script; this is
normal, and is just a side-effect of `numba` caching and pre-compiling some of the functions
used to drive the parcel model simulation.

> [!NOTE]
> We provide `pixi` environments for Linux, MacOS (both Intel and Apple Silicon) and
> Windows, but we have never tried to run the model on a Windows computer so your mileage
> may vary. Contact the authors if you have any questions and we can try to support your
> use case.

Installation
------------

To get started with using `pyrcel`, complete the following steps:

1. Set up a new Python environment; we recommend using [mambaforge](https://conda-forge.org/miniforge/):
  
``` shell
  $ mamba create -n pyrcel_quick_start python=3.11
```

2. Activate the new Python environment and install the model and its dependencies. If you install the published version from PyPi (_recommended_), then you also need to install [Assimulo](http://www.jmodelica.org/assimulo) using the Mamba package manager - but no other manual dependency installation is necessary:
  
``` shell
  $ mamba activate pyrcel_quick_start
  $ pip install pyrcel
  $ mamba install -c conda-forge assimulo
```

3. Run a test simulation using the CLI tool and a sample YAML file from **pyrcel/examples/\*.yml** (you may want to clone the repository or download them locally):
  
``` shell
  $ run_parcel simple.yml
```

* Visualize the output NetCDF (should be in the directory you ran the CLI tool, at **output/simple.nc**)

That's it! You should be able to import `pyrcel` into any script or program running in the
environment you created.


Requirements
------------

**Required**

* Python >= 3.8
* [numba](http://numba.pydata.org)
* [NumPy](http://www.numpy.org)
* [SciPy](http://www.scipy.org)
* [pandas](http://pandas.pydata.org)
* [xarray](http://xarray.pydata.org/en/stable/)
* [PyYAML](http://pyyaml.org/)

Additionally, the following packages are used for better numerics (ODE solving)

* [Assimulo](http://www.jmodelica.org/assimulo)

The easiest way to satisfy the basic requirements for building and running the
model is to use the [Anaconda](http://continuum.io/downloads) scientific Python
distribution. Alternatively, a
[miniconda environment](http://conda.pydata.org/docs/using/envs.html) is
provided to quickly set-up and get running the model. Assimulo's dependency on
the SUNDIALS library makes it a little bit tougher to install in an automated
fashion, so it has not been included in the automatic setup provided here; you
should refer to [Assimulo's documentation](http://www.jmodelica.org/assimulo_home/installation.html)
for more information on its installation process. Note that many components of
the model and package can be used without Assimulo.

Development
-----------

[http://github.com/darothen/pyrcel]()

Please fork this repository if you intend to develop the model further so that the
code's provenance can be maintained.

License / Usage
---------------

[All scientific code should be licensed](http://www.astrobetter.com/the-whys-and-hows-of-licensing-scientific-code/). This code is released under the New BSD (3-clause) [license](LICENSE.md).

You are free to use this code however you would like.
If you use this for any scientific work resulting in a publication or citation, please
cite our original publication detailing the model, and let the authors know:

```
@article { 
      author = "Daniel Rothenberg and Chien Wang",
      title = "Metamodeling of Droplet Activation for Global Climate Models",
      journal = "Journal of the Atmospheric Sciences",
      year = "2016",
      publisher = "American Meteorological Society",
      address = "Boston MA, USA",
      volume = "73",
      number = "3",
      doi = "10.1175/JAS-D-15-0223.1",
      pages= "1255 - 1272",
      url = "https://journals.ametsoc.org/view/journals/atsc/73/3/jas-d-15-0223.1.xml"
}
```


[author_email]: mailto:daniel@danielrothenberg.com
[nenes2001]: http://nenes.eas.gatech.edu/Preprints/KinLimitations_TellusPP.pdf
[pk2007]: http://www.atmos-chem-phys.net/7/1961/2007/acp-7-1961-2007.html
[hplgit]: https://github.com/hplgit/odespy
[diffrax]: https://docs.kidger.site/diffrax/
