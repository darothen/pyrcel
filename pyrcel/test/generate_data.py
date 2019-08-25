""" Generate test case data for future reference.
"""
import os
import pickle
from itertools import product

import numpy as np

from pyrcel import thermo

REFERENCE_FN = "results.dict"


def generate_reference(overwrite=False):

    results = dict()
    results["temperatures"] = np.linspace(233, 333, 10)  # K
    results["pressures"] = np.linspace(100, 1050, 10)  # hPa
    results["radii"] = np.logspace(-3, 1, 10)  # microns
    results["densities"] = np.linspace(0.5, 1.5, 10)
    results["dry_radii"] = np.logspace(-8, -6, 5)
    results["kappas"] = np.logspace(-2, 0, 5)
    results["r_over_r_dry"] = np.logspace(0.1, 3, 5)

    print("dv_cont", end=", ")
    results["dv_cont"] = [
        thermo.dv_cont(T, P * 100)
        for T, P in product(results["temperatures"], results["pressures"])
    ]
    print("(%d cases)" % len(results["dv_cont"]))

    print("dv", end=", ")
    results["dv"] = [
        thermo.dv(T, r * 1e-6, P * 100)
        for T, r, P in product(
            results["temperatures"], results["radii"], results["pressures"]
        )
    ]
    print("(%d cases)" % len(results["dv"]))

    print("rho_air", end=", ")
    results["rho_air"] = [
        thermo.rho_air(T, P * 100)
        for T, P in product(results["temperatures"], results["pressures"])
    ]
    print("(%d cases)" % len(results["rho_air"]))

    print("es", end=", ")
    results["es"] = [thermo.es(T - 273.15) for T in results["temperatures"]]
    print("(%d cases)" % len(results["es"]))

    print("ka_cont", end=", ")
    results["ka_cont"] = [thermo.ka_cont(T) for T in results["temperatures"]]
    print("(%d cases)" % len(results["ka_cont"]))

    print("ka", end=", ")
    results["ka"] = [
        thermo.ka(T, rho, r * 1e-6)
        for T, rho, r in product(
            results["temperatures"], results["densities"], results["radii"]
        )
    ]
    print("(%d cases)" % len(results["ka"]))

    print("sigma_w", end=", ")
    results["sigma_w"] = [thermo.sigma_w(T) for T in results["temperatures"]]
    print("(%d cases)" % len(results["sigma_w"]))

    print("Seq", end=", ")
    results["Seq_approx"] = [
        thermo.Seq(f * r_dry * 1e-6, r_dry * 1e-6, T, kappa)
        for f, r_dry, T, kappa in product(
            results["r_over_r_dry"],
            results["dry_radii"],
            results["temperatures"],
            results["kappas"],
        )
    ]
    results["Seq_exact"] = [
        thermo.Seq_approx(f * r_dry * 1e-6, r_dry * 1e-6, T, kappa)
        for f, r_dry, T, kappa in product(
            results["r_over_r_dry"],
            results["dry_radii"],
            results["temperatures"],
            results["kappas"],
        )
    ]
    print("(%d cases)" % (2 * len(results["Seq_exact"]),))

    if (not os.path.exists(REFERENCE_FN)) or overwrite:
        with open(REFERENCE_FN, "wb") as f:
            pickle.dump(results, f)


if __name__ == "__main__":

    generate_reference(True)
