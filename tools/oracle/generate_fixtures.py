"""Generate golden reference fixtures from the ``master`` parcel model.

This script is the *oracle* half of the v2 validation harness (design doc §7.0).
It must be run in an environment where the **current** (``master``) model works,
i.e. with ``numba`` and Assimulo/SUNDIALS ``CVode`` available. For this repo that
is the ``py311`` pixi environment:

    .pixi/envs/py311/bin/python tools/oracle/generate_fixtures.py

For each scenario in ``tests/scenarios.py`` it writes one ``.npz`` into
``tests/fixtures/oracle/`` containing:

  * aerosol arrays        : r_drys, kappas, Nis
  * equilibration fixture : y0  (initial state vector from ``_setup_run``)
  * RHS fixtures          : rhs_Y, rhs_dYdt evaluated by the numba ``parcel_ode_sys``
                            at physical (on-trajectory) + randomized states
  * trajectory fixture    : traj_t, traj_X (full CVode solution)
  * derived scalars       : smax, t_smax, z_smax
  * activation diagnostics: per-species eq/kn/alpha/phi at S_max

plus a ``manifest.json`` recording provenance (git commit, package versions,
solver tolerances) so the fixtures are reproducible and auditable.

The v2 (JAX/diffrax) test suite consumes these frozen fixtures and never needs
``numba`` or Assimulo at test time.
"""

from __future__ import annotations

import json
import subprocess
import sys
import warnings
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

warnings.filterwarnings("ignore")

REPO_ROOT = Path(__file__).resolve().parents[2]
FIXTURE_DIR = REPO_ROOT / "tests" / "fixtures" / "oracle"

# Import the shared scenario definitions directly (avoid package ambiguity).
sys.path.insert(0, str(REPO_ROOT / "tests"))
import scenarios as scn  # noqa: E402

# Number of sampled states for the RHS-equivalence fixtures.
N_PHYSICAL_STATES = 8
N_RANDOM_STATES = 16
RANDOM_SEED = 20260614


def _git_commit() -> str:
    try:
        return (
            subprocess.check_output(
                ["git", "rev-parse", "HEAD"], cwd=str(REPO_ROOT)
            )
            .decode()
            .strip()
        )
    except Exception:
        return "unknown"


def _versions() -> dict:
    info = {"python": sys.version.split()[0]}
    for mod in ("pyrcel", "numpy", "numba", "assimulo", "scipy"):
        try:
            m = __import__(mod)
            info[mod] = getattr(m, "__version__", "unknown")
        except Exception as exc:  # pragma: no cover - diagnostic only
            info[mod] = f"unavailable ({exc})"
    return info


def _sample_physical_states(X: np.ndarray, smax_idx: int, n: int) -> np.ndarray:
    """Pick ``n`` representative on-trajectory states, always including S_max."""
    m = X.shape[0]
    idx = np.unique(
        np.concatenate(
            [np.linspace(0, m - 1, n).round().astype(int), np.array([smax_idx])]
        )
    )
    return X[idx]


def _randomized_states(
    base_states: np.ndarray, r_drys: np.ndarray, n: int, rng: np.random.Generator
) -> np.ndarray:
    """Build randomized-but-valid states to stress the RHS off the trajectory."""
    n_state = 7  # N_STATE_VARS
    out = []
    attempts = 0
    while len(out) < n and attempts < 50 * n:
        attempts += 1
        base = base_states[rng.integers(len(base_states))].copy()
        y = base.copy()
        y[2] = base[2] + rng.uniform(-3.0, 3.0)  # T
        y[1] = base[1] * rng.uniform(0.97, 1.03)  # P
        y[3] = base[3] * rng.uniform(0.95, 1.05)  # wv
        y[6] = base[6] + rng.uniform(-5e-3, 5e-3)  # S
        rs = base[n_state:] * np.exp(rng.normal(0.0, 0.4, size=r_drys.shape))
        rs = np.maximum(rs, r_drys * 1.0001)  # keep r > r_dry (Seq well-defined)
        y[n_state:] = rs
        out.append(y)
    return np.asarray(out)


def generate_one(scenario: dict, rng: np.random.Generator) -> dict:
    import pyrcel as pm  # noqa: F401  (ensures package import)
    import pyrcel.constants as c
    from pyrcel._parcel_aux_numba import parcel_ode_sys
    from pyrcel.activation import binned_activation

    name = scenario["name"]
    ic = scenario["initial"]
    run_kw = scenario["run"]

    aerosols = scn.build_aerosols(scenario)
    model = pm.ParcelModel(
        aerosols,
        V=ic["V"],
        T0=ic["T0"],
        S0=ic["S0"],
        P0=ic["P0"],
        accom=ic["accom"],
        console=False,
    )

    nr = model._nr
    r_drys = np.asarray(model._r_drys, dtype=np.float64)
    kappas = np.asarray(model._kappas, dtype=np.float64)
    Nis = np.asarray(model._Nis, dtype=np.float64)
    y0 = np.asarray(model.y0, dtype=np.float64)

    parcel_df, aerosol_dfs = model.run(
        run_kw["t_end"],
        output_dt=run_kw["output_dt"],
        solver_dt=run_kw["solver_dt"],
        solver="cvode",
        output_fmt="dataframes",
        terminate=run_kw["terminate"],
        terminate_depth=run_kw["terminate_depth"],
        max_steps=run_kw["max_steps"],
    )

    traj_t = np.asarray(model.time, dtype=np.float64)
    traj_X = np.asarray(model.x, dtype=np.float64)

    s_idx = int(c.STATE_VAR_MAP["S"])
    z_idx = int(c.STATE_VAR_MAP["z"])
    smax_idx = int(np.argmax(traj_X[:, s_idx]))
    smax = float(traj_X[smax_idx, s_idx])
    t_smax = float(traj_t[smax_idx])
    z_smax = float(traj_X[smax_idx, z_idx])
    T_smax = float(traj_X[smax_idx, c.STATE_VAR_MAP["T"]])

    # --- RHS fixtures: evaluate the numba derivative at many states ---------------
    rhs_args = (nr, r_drys, Nis, float(ic["V"]), kappas, float(ic["accom"]))
    phys = _sample_physical_states(traj_X, smax_idx, N_PHYSICAL_STATES)
    rand = _randomized_states(phys, r_drys, N_RANDOM_STATES, rng)
    Y = np.vstack([phys, rand]).astype(np.float64)

    dYdt = np.empty_like(Y)
    for i in range(Y.shape[0]):
        dYdt[i] = np.asarray(parcel_ode_sys(Y[i], 0.0, *rhs_args), dtype=np.float64)
    if not np.all(np.isfinite(dYdt)):
        raise RuntimeError(f"non-finite RHS output for scenario {name!r}")

    # --- Activation diagnostics at S_max -----------------------------------------
    act_species, act_eq, act_kn, act_alpha, act_phi = [], [], [], [], []
    for aer in aerosols:
        rs_smax = np.asarray(aerosol_dfs[aer.species].iloc[smax_idx].values)
        eq, kn, alpha, phi = binned_activation(smax, T_smax, rs_smax, aer)
        act_species.append(aer.species)
        act_eq.append(float(eq))
        act_kn.append(float(kn))
        act_alpha.append(float(alpha))
        act_phi.append(float(phi))

    npz = dict(
        name=np.array(name),
        nr=np.array(nr),
        V=np.array(float(ic["V"])),
        accom=np.array(float(ic["accom"])),
        T0=np.array(float(ic["T0"])),
        S0=np.array(float(ic["S0"])),
        P0=np.array(float(ic["P0"])),
        r_drys=r_drys,
        kappas=kappas,
        Nis=Nis,
        y0=y0,
        rhs_Y=Y,
        rhs_dYdt=dYdt,
        n_physical=np.array(phys.shape[0]),
        n_random=np.array(rand.shape[0]),
        traj_t=traj_t,
        traj_X=traj_X,
        smax=np.array(smax),
        t_smax=np.array(t_smax),
        z_smax=np.array(z_smax),
        T_smax=np.array(T_smax),
        act_species=np.array(act_species),
        act_eq=np.array(act_eq, dtype=np.float64),
        act_kn=np.array(act_kn, dtype=np.float64),
        act_alpha=np.array(act_alpha, dtype=np.float64),
        act_phi=np.array(act_phi, dtype=np.float64),
    )

    summary = dict(
        name=name,
        nr=int(nr),
        n_steps=int(traj_X.shape[0]),
        smax=smax,
        t_smax=t_smax,
        z_smax=z_smax,
        activation={
            sp: {"eq": eq, "kn": kn}
            for sp, eq, kn in zip(act_species, act_eq, act_kn)
        },
    )
    return {"npz": npz, "summary": summary}


def main() -> int:
    import pyrcel.integrator as integ

    FIXTURE_DIR.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(RANDOM_SEED)

    summaries = []
    for scenario in scn.SCENARIOS:
        name = scenario["name"]
        print(f"[oracle] generating {name} ...", flush=True)
        result = generate_one(scenario, rng)
        out_path = FIXTURE_DIR / f"{name}.npz"
        np.savez_compressed(out_path, **result["npz"])
        summaries.append(result["summary"])
        s = result["summary"]
        print(
            f"         -> {out_path.name}: nr={s['nr']} steps={s['n_steps']} "
            f"Smax={s['smax']:.6e}",
            flush=True,
        )

    manifest = dict(
        generated_at=datetime.now(timezone.utc).isoformat(),
        git_commit=_git_commit(),
        versions=_versions(),
        solver=dict(
            backend="cvode (Assimulo/SUNDIALS, BDF)",
            rtol=integ.state_rtol,
            atol_state=integ.state_atol,
            atol_radius=1e-12,
        ),
        rhs=dict(
            n_physical=N_PHYSICAL_STATES,
            n_random=N_RANDOM_STATES,
            random_seed=RANDOM_SEED,
        ),
        scenarios=summaries,
    )
    with open(FIXTURE_DIR / "manifest.json", "w") as f:
        json.dump(manifest, f, indent=2)
    print(f"[oracle] wrote {len(summaries)} fixtures + manifest to {FIXTURE_DIR}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
