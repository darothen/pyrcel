"""pyrcel: 0-D adiabatic cloud parcel model (v2, JAX/diffrax backend).

Public API
----------
The primary entry points are:

* :class:`~pyrcel.model.ParcelModelJAX` — the v2 parcel model
* :class:`~pyrcel.aerosol.AerosolSpecies` — aerosol population container
* :class:`~pyrcel.distributions.Lognorm` — log-normal size distribution
* :class:`~pyrcel.updraft.ConstantV`, :class:`~pyrcel.updraft.InterpolatedUpdraft`

Legacy (NumPy/numba) implementations are preserved in :mod:`pyrcel.legacy`
for reference and cross-checking — not for production use.
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
    "ParcelModelJAX": "pyrcel.model",
    "run_updraft_ensemble": "pyrcel.ensemble",
    "smax_nact_ensemble": "pyrcel.ensemble",
    "sample_gaussian_updrafts": "pyrcel.ensemble",
    # Legacy class kept for the run_parcel CLI and any callers using the old API.
    "ParcelModel": "pyrcel.legacy.parcel",
    "run_model": "pyrcel.legacy.driver",
    "iterate_runs": "pyrcel.legacy.driver",
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
