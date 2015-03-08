# Parcel Model

![sample parcel model run](doc/figs/model_example.png)

This is an implementation of a simple, adiabatic cloud parcel model for use in aerosol-cloud interaction studies. It is based on the model used by [Nenes et al (2001)][Nenes2001], but with several key modifications:

* Implementation of kappa-Kohler theory for condensation physics ([Petters and Kreidenweis, 2007)][pk2007]
* Extension of model to handle arbitrary sectional representations of aerosol populations, based on user-controlled empirical or parameterized size distributions
* Improved, modular numerical framework for integrating the model, including bindings to several different stiff integrators:
 * `lsoda` - [scipy ODEINT wrapper](http://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.odeint.html)
 * `vode, lsode*, lsoda*` - ODEPACK via [odespy][hplgit]
 * `cvode` - SUNDIALS via [Assimulo](http://www.jmodelica.org/assimulo_home/index.html#)

among other details. It also includes a library of droplet activation routines and scripts/notebooks for evaluating those schemes against equivalent calculations done with the parcel model.

Updated code can be found the project [github repository](https://github.com/darothen/parcel_model). If you'd like to use this code or have any questions about it, please [contact the author][author_email]. In particular, if you use this code for research purposes, be sure to carefully read through the model and ensure that you have tweaked/configured it for your purposes (i.e., modifying the accomodation coefficient); other derived quantities). 

[Detailed documentation is available](http://mit.edu/~darothen/parcel_model/), including a [scientific description](http://mit.edu/~darothen/parcel_model/sci_descr.html), [installation details](http://mit.edu/~darothen/parcel_model/install.html), and a [basic example](http://mit.edu/~darothen/parcel_model/examples/basic_run.html) which produces a figure like the plot at the top of this page. 

## Development

[http://github.com/darothen/parcel_model]()

Please fork this repository if you intend to develop the model further so that the code's provenance can be maintained.

## License

[All scientific code should be licensed](http://www.astrobetter.com/the-whys-and-hows-of-licensing-scientific-code/). This code is released under the New BSD (3-clause) [license](LICENSE.md).

[author_email]: mailto:darothen@mit.edu

[nenes2001]: http://nenes.eas.gatech.edu/Preprints/KinLimitations_TellusPP.pdf
[pk2007]: http://www.atmos-chem-phys.net/7/1961/2007/acp-7-1961-2007.html
[hplgit]: https://github.com/hplgit/odespy
