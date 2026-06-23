"""Smoke-tests for the pyrcel top-level public API surface (PR #38).

Verifies that all documented names are importable from the top-level package,
that legacy-only names are NOT in the top-level namespace, and that lazy
attributes resolve to the correct objects.
"""

from __future__ import annotations

import pytest

# --- v2 names that must be importable from `pyrcel` directly --------------------


def test_parcel_model_is_v2():
    import pyrcel
    from pyrcel.model import ParcelModel

    assert pyrcel.ParcelModel is ParcelModel


def test_parcel_model_jax_alias():
    import pyrcel
    from pyrcel.model import ParcelModel

    assert pyrcel.ParcelModelJAX is ParcelModel


def test_aerosol_species():
    from pyrcel import AerosolSpecies
    from pyrcel.aerosol import AerosolSpecies as _Src

    assert AerosolSpecies is _Src


def test_distributions():
    from pyrcel import Lognorm, MultiModeLognorm
    from pyrcel.distributions import Lognorm as _L
    from pyrcel.distributions import MultiModeLognorm as _ML

    assert Lognorm is _L
    assert MultiModeLognorm is _ML


def test_updraft_classes():
    from pyrcel import AbstractUpdraft, ConstantV, InterpolatedUpdraft, as_updraft
    from pyrcel.updraft import (
        AbstractUpdraft as _AU,
    )
    from pyrcel.updraft import (
        ConstantV as _CV,
    )
    from pyrcel.updraft import (
        InterpolatedUpdraft as _IU,
    )
    from pyrcel.updraft import (
        as_updraft as _au,
    )

    assert ConstantV is _CV
    assert InterpolatedUpdraft is _IU
    assert AbstractUpdraft is _AU
    assert as_updraft is _au


def test_activation_exports():
    import pyrcel
    from pyrcel.activation import ARG2000 as _ARG
    from pyrcel.activation import ActivationScheme as _AS

    assert pyrcel.ARG2000 is _ARG
    assert pyrcel.ActivationScheme is _AS


def test_ensemble_exports():
    import pyrcel
    from pyrcel.ensemble import (
        run_updraft_ensemble as _rue,
    )
    from pyrcel.ensemble import (
        sample_gaussian_updrafts as _sgu,
    )
    from pyrcel.ensemble import (
        smax_nact_ensemble as _sne,
    )

    assert pyrcel.run_updraft_ensemble is _rue
    assert pyrcel.smax_nact_ensemble is _sne
    assert pyrcel.sample_gaussian_updrafts is _sgu


# --- legacy names must NOT be at the top level ----------------------------------


def test_run_model_not_in_top_level():
    import pyrcel

    with pytest.raises(AttributeError):
        _ = pyrcel.run_model


def test_iterate_runs_not_in_top_level():
    import pyrcel

    with pytest.raises(AttributeError):
        _ = pyrcel.iterate_runs


def test_legacy_numerical_oracles_accessible():
    """thermo and activation remain importable as cross-check references."""
    import pyrcel.legacy.activation
    import pyrcel.legacy.thermo

    assert hasattr(pyrcel.legacy.thermo, "Seq")
    assert hasattr(pyrcel.legacy.activation, "arg2000")
