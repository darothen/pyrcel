# Oracle fixture generation

This directory holds the **oracle** half of the v2 validation harness described in
[`docs/design/jax-diffrax-migration.md`](../../docs/design/jax-diffrax-migration.md) (§7).

The goal: freeze golden reference data from the **current `master`** model
(numba + Assimulo/SUNDIALS `CVode`) so the v2 (JAX + `diffrax`) model can be validated
against it **without** numba or Assimulo ever being needed at test time.

## What gets generated

`generate_fixtures.py` runs every scenario in [`tests/scenarios.py`](../../tests/scenarios.py)
and writes, per scenario, a compressed `.npz` into `tests/fixtures/oracle/` containing:

- aerosol arrays (`r_drys`, `kappas`, `Nis`) and the equilibrated initial state `y0`;
- **RHS fixtures** (`rhs_Y`, `rhs_dYdt`): the numba `parcel_ode_sys` evaluated at a set of
  on-trajectory ("physical") and randomized states — the cornerstone of the RHS-equivalence
  test (design doc §7.2);
- the full CVode **trajectory** (`traj_t`, `traj_X`);
- derived scalars (`smax`, `t_smax`, `z_smax`, `T_smax`);
- per-species **activation diagnostics** at `S_max` (`act_eq`, `act_kn`, `act_alpha`, `act_phi`).

It also writes `manifest.json` recording provenance: git commit, package versions, and the
exact solver tolerances used.

## How to regenerate

Requires an environment with the `master` model working (numba + Assimulo). In this repo
that is the `py311` pixi environment:

```shell
# from the repo root
<path-to>/.pixi/envs/py311/bin/python tools/oracle/generate_fixtures.py
```

(or `pixi run -e py311 python tools/oracle/generate_fixtures.py`)

Regenerate only when the oracle definition changes (new scenarios, or an intentional
change to the reference model). The fixtures are committed so day-to-day v2 development and
CI need neither numba nor SUNDIALS.

## Environments at a glance

| Purpose | Environment | Key deps |
| --- | --- | --- |
| **Oracle** (bake fixtures) | `py311` pixi env (or any `master` env) | numba, assimulo/CVode |
| **v2 dev/test** | JAX venv (e.g. `pip install -e ".[jax,test]"`) | jax, diffrax (→ optimistix, equinox, lineax) |

The v2 test suite (`tests/`) only needs the JAX venv; it reads the frozen fixtures.
