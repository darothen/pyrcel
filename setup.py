#! /usr/bin/env python
#
# Copyright (C) 2015 Daniel Rothenberg <darothen@mit.edu>

DESCRIPTION = "parcel_model: 0D adiabatic cloud parcel model"
LONG_DESCRIPTION = """\
This code implements a relatively simple, adiabatic cloud parcel model for studying aerosol
activation. It allows the user to peform and iterate parcel model experiments in a simple 
fashion through the use of an object-oriented implementation. Furthermore, interfaces for
several numerical solvers are built into the model so that users can choose whatever scheme
they would like to run with the model.
"""

DISTNAME = "parcel_model"
MAINTAINER = "Daniel Rothenberg"
AUTHOR = MAINTAINER
MAINTAINER_EMAIL = "darothen@mit.edu"
AUTHOR_EMAIL = MAINTAINER_EMAIL
URL = "https://github.com/darothen/parcel_model"
LICENSE = "New BSD (3-clause)"
DOWNLOAD_URL = "https://github.com/darothen/parcel_model"

execfile("parcel_model/version.py")
VERSION = __version__

try:
    from setuptools import setup
    HAS_SETUPTOOLS = True
except ImportError:
    from distutils.core import setup

# Must have Cython!!
from Cython.Distutils import build_ext
from distutils.extension import Extension
import numpy

extensions = [
    Extension("parcel_model.parcel_aux", 
              ["parcel_model/parcel_aux.pyx"],
              include_dirs=[numpy.get_include(), ],
              #extra_compile_args=['-fopenmp', ],
              #extra_link_args=['-fopenmp', ],
            ),
   #Extension("parcel_aux_bin", ["parcel_aux_bin.pyx"],
   #          libraries=cgsl.get_libraries(),
   #          library_dirs=[cgsl.get_library_dir(), ],
   #          include_dirs=[numpy.get_include(), cgsl.get_cython_include_dir()],
   #          ),
]

if __name__ == "__main__":

    setup(
        name = DISTNAME,
        author = AUTHOR,
        author_email = AUTHOR_EMAIL,
        maintainer = MAINTAINER,
        maintainer_email = MAINTAINER_EMAIL,
        description = DESCRIPTION,
        long_description = LONG_DESCRIPTION,
        license = LICENSE,
        url = URL,
        version = VERSION,
        download_url = DOWNLOAD_URL,

        install_requires = ['numpy', 'scipy', 'pandas', ],
        packages = ["parcel_model", ],
        ext_modules = extensions,
        cmdclass = {'build_ext': build_ext},

        classifiers=[
            'Development Status :: 5 - Production/Stable', 
            'Environment :: Console',
            'Intended Audience :: Science/Research',
            'License :: OSI Approved :: MIT License',
            'Natural Language :: English',
            'Operating System :: Unix',
            'Programming Language :: Python :: 2.6',
            'Programming Language :: Python :: 2.7',
            'Programming Language :: Fortran',
            'Topic :: Scientific/Engineering :: Atmospheric Science',
        ],
    )
