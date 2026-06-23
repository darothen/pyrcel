"""Benchmark the v2 (JAX + diffrax) parcel-model core.

Run with the **v2 environment** (JAX + diffrax; no numba/Assimulo needed):

    .venv/bin/python tools/bench/bench_v2.py [scenario ...] [--batch N] [--reps R]

For each scenario it times:

* **equilibration** (``optimistix`` root-finds, jitted),
* **integration** -- cold (first call, includes XLA compile) and warm (median), and
* optionally a **vmap ensemble** over ``N`` perturbed initial conditions, reporting
  amortized per-solve cost (the batching/GPU advantage).

Initial state ``y0`` and the output time grid are taken from the frozen oracle
fixtures (``tests/fixtures/oracle``) so the numbers line up 1:1 with
``bench_master.py`` (identical horizon, identical tolerances ``rtol=1e-7`` + the
per-component ``atol`` vector). See ``tools/bench/README.md``.
"""

from __future__ import annotations

import argparse
import statistics
import sys
import time
import warnings
from pathlib import Path

warnings.filterwarnings("ignore")

import numpy as np

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "tests"))
FIXTURE_DIR = REPO / "tests" / "fixtures" / "oracle"

import jax  # noqa: E402

jax.config.update("jax_enable_x64", True)

import jax.numpy as jnp  # noqa: E402

from pyrcel.equilibrate_jax import equilibrate_initial_state  # noqa: E402
from pyrcel.integrator_diffrax import _solve, atol_vector, integrate_parcel  # noqa: E402
from pyrcel.parcel_aux_jax import N_STATE_VARS  # noqa: E402


def _median(fn, reps: int) -> float:
    samples = []
    for _ in range(reps):
        t0 = time.perf_counter()
        fn()
        samples.append(time.perf_counter() - t0)
    return statistics.median(samples)


def bench(name: str, batch: int, reps: int) -> None:
    d = np.load(FIXTURE_DIR / f"{name}.npz")
    y0 = jnp.asarray(d["y0"])
    ts = jnp.asarray(d["traj_t"])
    nr = int(y0.shape[0] - N_STATE_VARS)
    args = (
        jnp.asarray(d["r_drys"]),
        jnp.asarray(d["Nis"]),
        jnp.asarray(d["kappas"]),
        float(d["accom"]),
        float(d["V"]),
    )

    # Integration: cold (incl. compile) then warm median.
    t0 = time.perf_counter()
    integrate_parcel(y0, args, ts).ys.block_until_ready()
    cold = time.perf_counter() - t0
    warm = _median(lambda: integrate_parcel(y0, args, ts).ys.block_until_ready(), reps)

    # Equilibration (jitted; warm).
    eq = jax.jit(
        lambda: equilibrate_initial_state(
            float(d["T0"]),
            float(d["S0"]),
            float(d["P0"]),
            jnp.asarray(d["r_drys"]),
            jnp.asarray(d["kappas"]),
            jnp.asarray(d["Nis"]),
        )
    )
    eq().block_until_ready()
    eq_warm = _median(lambda: eq().block_until_ready(), reps)

    line = (
        f"{name:16s} nr={nr:3d} npts={len(ts):4d} | "
        f"equil={eq_warm * 1e3:6.2f}ms | "
        f"integ cold={cold * 1e3:8.1f}ms warm={warm * 1e3:8.2f}ms"
    )

    if batch > 0:
        atol = atol_vector(nr)
        rng = np.random.default_rng(0)
        Y0 = jnp.asarray(
            np.asarray(d["y0"])[None, :] * (1.0 + 1e-4 * rng.standard_normal((batch, y0.shape[0])))
        )
        batched = jax.jit(
            lambda Y: jax.vmap(lambda yy: _solve(yy, args, ts, 1e-7, atol, None, 100_000).ys)(Y)
        )
        t0 = time.perf_counter()
        batched(Y0).block_until_ready()
        b_compile = time.perf_counter() - t0
        b_warm = _median(lambda: batched(Y0).block_until_ready(), max(3, reps // 4))
        line += (
            f" | vmap{batch} compile={b_compile:5.1f}s "
            f"warm={b_warm * 1e3:8.1f}ms ({b_warm / batch * 1e3:6.3f}ms/solve)"
        )

    print(line, flush=True)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("scenarios", nargs="*", help="scenario names (default: all)")
    parser.add_argument("--batch", type=int, default=0, help="vmap ensemble size (0 = skip)")
    parser.add_argument("--reps", type=int, default=20, help="warm-timing repetitions")
    a = parser.parse_args()

    import scenarios as scn

    names = a.scenarios or [s["name"] for s in scn.SCENARIOS]
    print(f"[bench-v2] jax {jax.__version__} | device {jax.devices()[0]}", flush=True)
    for name in names:
        bench(name, a.batch, a.reps)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
