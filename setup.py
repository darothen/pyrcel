#! /usr/bin/env python

try:
    from setuptools import setup

    HAS_SETUPTOOLS = True
except ImportError:
    from distutils.core import setup

import os
import warnings
from textwrap import dedent

MAJOR, MINOR, MICRO = 1, 3, 1
DEV_ITER = 0
DEV = True

if DEV:
    VERSION = "{}.{}.dev{}".format(MAJOR, MINOR, DEV_ITER)
else:
    VERSION = "{}.{}.{}".format(MAJOR, MINOR, MICRO)

extensions = [
    # Numba AOT extension module
    # auxcc.distutils_extension(),
]


def _write_version_file():
    fn = os.path.join(os.path.dirname(__file__), "pyrcel", "version.py")

    version_str = dedent(
        """
        __version__ = '{}'
        """
    )

    # Write version file
    with open(fn, "w") as version_file:
        version_file.write(version_str.format(VERSION))


# Write version and install
_write_version_file()

setup(
    name="pyrcel",
    author="Daniel Rothenberg",
    author_email="daniel@danielrothenberg.com",
    maintainer="Daniel Rothenberg",
    maintainer_email="daniel@danielrothenberg.com",
    description="pyrcel: 0D adiabatic cloud parcel model",
    long_description="""
        This code implements a relatively simple, adiabatic cloud parcel model for studying aerosol
        activation. It allows the user to peform and iterate parcel model experiments in a simple
        fashion through the use of an object-oriented implementation. Furthermore, interfaces for
        several numerical solvers are built into the model so that users can choose whatever scheme
        they would like to run with the model.
    """,
    license="New BSD (3-clause)",
    url="https://github.com/darothen/pyrcel",
    version=VERSION,
    download_url="https://github.com/darothen/pyrcel",
    # TODO: Update install requirements and corresponding documentation
    install_requires=[
        "numba",
        "numpy",
        "pandas",
        "pyyaml",
        "scipy",
        "setuptools",
        "xarray",
    ],
    packages=["pyrcel"],
    package_data={"pyrcel": ["data/std_atm.csv"]},
    scripts=["scripts/run_parcel"],
    ext_modules=extensions,
    classifiers=[
        "Development Status :: 5 - Production/Stable",
        "Environment :: Console",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: BSD License",
        "Natural Language :: English",
        "Operating System :: Unix",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Topic :: Scientific/Engineering :: Atmospheric Science",
    ],
)
