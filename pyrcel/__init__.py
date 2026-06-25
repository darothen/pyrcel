"""pyrcel: 0-D adiabatic cloud parcel model (v2, JAX/diffrax backend).

Public API
----------
The stable v2 entry points are:

* [ParcelModel][pyrcel.model.ParcelModel] — v2 parcel model (JAX/diffrax, differentiable)
* [AerosolSpecies][pyrcel.aerosol.AerosolSpecies] — aerosol population container
* [Lognorm][pyrcel.distributions.Lognorm], [MultiModeLognorm][pyrcel.distributions.MultiModeLognorm]
* [ConstantV][pyrcel.updraft.ConstantV], [InterpolatedUpdraft][pyrcel.updraft.InterpolatedUpdraft],
  [as_updraft][pyrcel.updraft.as_updraft]
* [ModelOutput][pyrcel.model_output.ModelOutput] — structured output returned by
  ``ParcelModel.run(mode='full')``; exposes ``.to_pandas()``, ``.to_polars()``,
  ``.to_xarray()``, ``.to_netcdf()``, ``.to_csv()``, ``.to_parquet()``
* [ARG2000][pyrcel.activation.ARG2000], [MBN2014][pyrcel.activation.MBN2014],
  [ActivationScheme][pyrcel.activation.ActivationScheme]
* [run_updraft_ensemble][pyrcel.ensemble.run_updraft_ensemble],
  [smax_nact_ensemble][pyrcel.ensemble.smax_nact_ensemble]

The legacy numerical oracles (NumPy/SciPy) are preserved in `pyrcel.legacy` for
cross-checking only — `pyrcel.legacy.thermo` and `pyrcel.legacy.activation`
are the primary reference implementations. The legacy parcel model (CVode/Assimulo)
has been removed; a v2 CLI is tracked in issue #60.
"""

from importlib.metadata import version as _version

try:
    __version__ = _version("pyrcel")
except Exception:
    __version__ = "local"

__author__ = "Daniel Rothenberg <daniel@danielrothenberg.com>"

from .aerosol import AerosolSpecies
from .distributions import Lognorm, MultiModeLognorm
from .updraft import AbstractUpdraft, ConstantV, InterpolatedUpdraft, as_updraft

_LAZY_ATTRS = {
    # v2 primary
    "ParcelModel": "pyrcel.model",
    "ModelOutput": "pyrcel.model_output",
    # backward-compat alias for ParcelModel
    "ParcelModelJAX": "pyrcel.model",
    # activation (JAX-heavy; imported lazily)
    "ARG2000": "pyrcel.activation",
    "MBN2014": "pyrcel.activation",
    "ActivationScheme": "pyrcel.activation",
    # ensemble utilities
    "run_updraft_ensemble": "pyrcel.ensemble",
    "smax_nact_ensemble": "pyrcel.ensemble",
    "sample_gaussian_updrafts": "pyrcel.ensemble",
}


def __getattr__(name):
    module_path = _LAZY_ATTRS.get(name)
    if module_path is not None:
        import importlib

        module = importlib.import_module(module_path)
        return getattr(module, name)
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")


def __dir__():
    return sorted(list(globals()) + list(_LAZY_ATTRS))
