""" Test cases for thermodynamics module.

Most of these test cases just compare the current version of the code's results
for parameter sets versus a reference to serve as basic regression testing.

"""
from __future__ import absolute_import

import os, pickle
import unittest

from itertools import product

import numpy as np
from numpy.testing import assert_allclose

from .. thermo import *
from . generate_data import REFERENCE_FN

class TestThermoTestCases(unittest.TestCase):

    def setUp(self):

        with open(REFERENCE_FN, 'rb') as f:
            self.reference = pickle.load(f)

        self.temperatures = self.reference['temperatures']
        self.pressures = self.reference['pressures']
        self.radii = self.reference['radii']

    def test_dv_cont(self):
        result = [ dv_cont(T, P*100) for T, P in \
                    product(self.temperatures, self.pressures) ]

        assert_allclose(self.reference['dv_cont'], result)

    def test_dv(self):
        result = [ dv(T, r*1e-6, P*100) for T, r, P in \
                    product(self.temperatures, self.radii, self.pressures) ]

        assert_allclose(self.reference['dv'], result)
