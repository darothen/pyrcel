"""Benchmark the master (numba + Assimulo/CVode) parcel model.

Run with the **oracle environment** that has numba + Assimulo/SUNDIALS available
(for this repo, the ``py311`` pixi env):

    .pixi/envs/py311/bin/python tools/bench/bench_master.py [scenario ...] [--reps R]

For each scenario it times:

* **equilibration** -- the ``ParcelModel`` constructor (``scipy.bisect``/``fminbound``),
* **integration** -- ``model.run(solver="cvode")``, first call (includes the numba
  JIT compile) and warm (median).

Same scenarios, horizon, and solver tolerances as ``bench_v2.py`` so the two are
directly comparable. See ``tools/bench/README.md``.
"""

from __future__ import annotations

import argparse
import statistics
import sys
import time
import warnings
from pathlib import Path

warnings.filterwarnings("ignore")

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "tests"))

import scenarios as scn  # noqa: E402

import pyrcel as pm  # noqa: E402


def _median(fn, reps: int) -> float:
    samples = []
    for _ in range(reps):
        t0 = time.perf_counter()
        fn()
        samples.append(time.perf_counter() - t0)
    return statistics.median(samples)


def bench(name: str, reps: int) -> None:
    sc = scn.get_scenario(name)
    ic = sc["initial"]
    run = sc["run"]
    aerosols = scn.build_aerosols(sc)

    def construct():
        return pm.ParcelModel(
            aerosols,
            V=ic["V"],
            T0=ic["T0"],
            S0=ic["S0"],
            P0=ic["P0"],
            accom=ic["accom"],
            console=False,
        )

    eq_warm = _median(construct, max(5, reps // 2))
    model = construct()

    run_kw = dict(
        output_dt=run["output_dt"],
        solver_dt=run["solver_dt"],
        solver="cvode",
        output_fmt="arrays",
        terminate=run["terminate"],
        terminate_depth=run["terminate_depth"],
        max_steps=run["max_steps"],
    )

    t0 = time.perf_counter()
    out = model.run(run["t_end"], **run_kw)
    first = time.perf_counter() - t0
    warm = _median(lambda: model.run(run["t_end"], **run_kw), reps)

    npts = out[0].shape[0]
    print(
        f"{name:16s} nr={model._nr:3d} npts={npts:4d} | "
        f"equil={eq_warm * 1e3:6.2f}ms | "
        f"run first={first * 1e3:8.1f}ms warm={warm * 1e3:8.2f}ms",
        flush=True,
    )


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("scenarios", nargs="*", help="scenario names (default: all)")
    parser.add_argument("--reps", type=int, default=8, help="warm-timing repetitions")
    a = parser.parse_args()

    names = a.scenarios or [s["name"] for s in scn.SCENARIOS]
    print(f"[bench-master] pyrcel {pm.__version__}", flush=True)
    for name in names:
        bench(name, a.reps)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
