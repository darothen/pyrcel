"""pyrcel: 0-D adiabatic cloud parcel model (v2, JAX/diffrax backend).

Public API
----------
The stable v2 entry points are:

* :class:`~pyrcel.model.ParcelModel` — v2 parcel model (JAX/diffrax, differentiable)
* :class:`~pyrcel.aerosol.AerosolSpecies` — aerosol population container
* :class:`~pyrcel.distributions.Lognorm`, :class:`~pyrcel.distributions.MultiModeLognorm`
* :class:`~pyrcel.updraft.ConstantV`, :class:`~pyrcel.updraft.InterpolatedUpdraft`,
  :func:`~pyrcel.updraft.as_updraft`
* :class:`~pyrcel.activation.ARG2000`, :class:`~pyrcel.activation.ActivationScheme`
* :func:`~pyrcel.ensemble.run_updraft_ensemble`,
  :func:`~pyrcel.ensemble.smax_nact_ensemble`

Legacy (NumPy/numba) implementations are preserved in :mod:`pyrcel.legacy` as a
cross-check oracle and for callers that have not yet migrated. They are no longer
exported from this top-level namespace:

* :class:`pyrcel.legacy.parcel.ParcelModel` — legacy CVode/Assimulo model
* :func:`pyrcel.legacy.driver.run_model`, :func:`pyrcel.legacy.driver.iterate_runs`
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
    # backward-compat alias for ParcelModel
    "ParcelModelJAX": "pyrcel.model",
    # activation (JAX-heavy; imported lazily)
    "ARG2000": "pyrcel.activation",
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
