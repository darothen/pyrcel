#! /usr/bin/env python

try:
    from setuptools import setup
    HAS_SETUPTOOLS = True
except ImportError:
    from distutils.core import setup

# Must have Cython!!
from Cython.Distutils import build_ext
from distutils.extension import Extension
import numpy

import os
import warnings
from textwrap import dedent

MAJOR, MINOR, MICRO = 1, 3, 1
DEV = False
VERSION = "{}.{}.{}".format(MAJOR, MINOR, MICRO)

# Correct versioning with git info if DEV
if DEV:
    import subprocess

    pipe = subprocess.Popen(
        ['git', "describe", "--always", "--match", "'v[0-9]*'"],
        stdout=subprocess.PIPE
    )
    so, err = pipe.communicate()

    if pipe.returncode != 0:
        # no git or something wrong with git (not in dir?)
        warnings.warn("WARNING: Couldn't identify git revision, using generic version string")
        VERSION += ".dev"
    else:
        git_rev = so.strip()
        git_rev = git_rev.decode('ascii')  # necessary for Python >= 3

        VERSION += ".dev-" + format(git_rev)

extensions = [
    Extension("pyrcel.parcel_aux",
              ["pyrcel/parcel_aux.pyx"],
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

def _write_version_file():

    fn = os.path.join(os.path.dirname(__file__), 'pyrcel', 'version.py')

    version_str = dedent("""
        __version__ = '{}'
        """)

    # Write version file
    with open(fn, 'w') as version_file:
        version_file.write(version_str.format(VERSION))

# Write version and install
_write_version_file()

setup(
    name = "pyrcel",
    author = "Daniel Rothenberg",
    author_email = "darothen@mit.edu",
    maintainer = "Daniel Rothenberg",
    maintainer_email = "darothen@mit.edu",
    description = "pyrcel: 0D adiabatic cloud parcel model",
    long_description = """
        This code implements a relatively simple, adiabatic cloud parcel model for studying aerosol
        activation. It allows the user to peform and iterate parcel model experiments in a simple
        fashion through the use of an object-oriented implementation. Furthermore, interfaces for
        several numerical solvers are built into the model so that users can choose whatever scheme
        they would like to run with the model.
    """,
    license = "New BSD (3-clause)",
    url = "https://github.com/darothen/pyrcel",
    version = VERSION,
    download_url = "https://github.com/darothen/pyrcel",

    # TODO: Update install requirements and corresponding documentation
    install_requires = ['numpy', 'scipy', 'pandas', 'future'],
    packages = ["pyrcel", ],
    package_data = {
        'pyrcel': ['data/std_atm.csv', ],
    },
    scripts = [
        'scripts/run_parcel',
    ],
    ext_modules = extensions,
    cmdclass = {'build_ext': build_ext},

    classifiers=[
        'Development Status :: 5 - Production/Stable',
        'Environment :: Console',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License',
        'Natural Language :: English',
        'Operating System :: Unix',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Fortran',
        'Topic :: Scientific/Engineering :: Atmospheric Science',
    ],
)
