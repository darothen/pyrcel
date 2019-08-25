""" Test cases for thermodynamics module.

Most of these test cases just compare the current version of the code's results
for parameter sets versus a reference to serve as basic regression testing.

"""
import pickle
import unittest
from itertools import product

from numpy.testing import assert_allclose

from .generate_data import REFERENCE_FN
from ..thermo import *


class TestThermoTestCases(unittest.TestCase):
    def setUp(self):

        with open(REFERENCE_FN, "rb") as f:
            self.reference = pickle.load(f)

        self.temperatures = self.reference["temperatures"]
        self.pressures = self.reference["pressures"]
        self.radii = self.reference["radii"]
        self.densities = self.reference["densities"]

    def test_dv_cont(self):
        result = [
            dv_cont(T, P * 100)
            for T, P in product(self.temperatures, self.pressures)
        ]

        assert_allclose(self.reference["dv_cont"], result)

    def test_dv(self):
        result = [
            dv(T, r * 1e-6, P * 100)
            for T, r, P in product(
                self.temperatures, self.radii, self.pressures
            )
        ]

        assert_allclose(self.reference["dv"], result)

    def test_rho_air(self):
        result = [
            rho_air(T, P * 100)
            for T, P in product(self.temperatures, self.pressures)
        ]

        assert_allclose(self.reference["rho_air"], result)

    def test_es(self):
        result = [es(T - 273.15) for T in self.temperatures]

        assert_allclose(self.reference["es"], result)

    def test_ka_cont(self):
        result = [ka_cont(T) for T in self.temperatures]

        assert_allclose(self.reference["ka_cont"], result)

    def test_ka(self):
        result = [
            ka(T, rho, r * 1e-6)
            for T, rho, r in product(
                self.temperatures, self.densities, self.radii
            )
        ]

        assert_allclose(self.reference["ka"], result)

    def sigma_w(self):
        result = [sigma_w(T) for T in self.temperatures]

        assert_allclose(self.reference["sigma_w"], result)
