# Parcel-model benchmarks (v2 diffrax vs. master numba/CVode)

Two scripts that time the new JAX/diffrax model against the original numba +
Assimulo/CVode model on the shared scenario matrix (`tests/scenarios.py`), so the
performance trade-off of the migration is documented honestly (design doc §7.5).

The comparison is **apples-to-apples**: identical scenarios, identical integration
horizon and output grid (taken from the frozen oracle fixtures), and identical solver
tolerances (`rtol=1e-7` and the per-component `atol` vector). The two models agree on
`S_max` to ~1e-5 and on activated fraction exactly (see `tests/test_trajectory.py`).

## Running

v2 (JAX + diffrax) -- in the v2 environment:

```bash
.venv/bin/python tools/bench/bench_v2.py                 # all scenarios
.venv/bin/python tools/bench/bench_v2.py simple_sulfate two_mode --batch 128
```

master (numba + Assimulo/CVode) -- in the oracle environment:

```bash
.pixi/envs/py311/bin/python tools/bench/bench_master.py
```

## What gets reported

- **equil**: initial-condition equilibration (v2: `optimistix` + `vmap`; master: the
  `ParcelModel` constructor's `scipy.bisect`/`fminbound`).
- **integ/run cold/first**: first call, including one-time compilation (XLA for v2,
  numba for master).
- **integ/run warm**: median wall-clock once compiled.
- **vmap{N}** (v2 only, with `--batch N`): one compiled solve over a batch of `N`
  perturbed parcels, plus amortized ms/solve -- the ensemble/GPU advantage that the
  numba/CVode stack cannot offer.

## Reference numbers (Apple Silicon, CPU, float64; indicative, not a guarantee)

| scenario | nr | master run (warm) | v2 integ (warm) | master equil | v2 equil |
| --- | --- | --- | --- | --- | --- |
| simple_sulfate | 50 | ~285 ms | ~166 ms | ~14 ms | ~0.5 ms |
| two_mode | 150 | ~2271 ms | ~774 ms | ~45 ms | ~1.3 ms |

One-time compile: master ~1.7-2.4 s (numba), v2 ~4.3 s (XLA, per distinct `nr`).
`vmap` ensemble (simple_sulfate, batch 128, CPU): ~37 ms/solve (~4.5x the single-run
throughput) after a ~12 s batched compile.

**Takeaways:** v2 is faster warm on a single CPU run and scales better with bin
count; equilibration is ~30x faster; batching (and GPU) multiplies throughput
further. The original wins only on a single one-shot run, due to its lower compile
latency. And the migration unlocks exact autodiff, which the numba/CVode stack
cannot provide at all.
```
