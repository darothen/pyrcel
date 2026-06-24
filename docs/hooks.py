"""mkdocs build hook: execute example scripts and capture output.

Set the environment variable ``DOCS_EXECUTE=1`` to run the scripts during
``mkdocs build``. Without it the hook only ensures the asset directories exist
(useful for ``mkdocs serve`` during local development, where pre-committed
placeholder assets are used instead).
"""

from __future__ import annotations

import os
import subprocess
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).parent.parent
DOCS_DIR = Path(__file__).parent
ASSETS_FIGURES = DOCS_DIR / "assets" / "figures"
ASSETS_OUTPUT = DOCS_DIR / "assets" / "output"

SCRIPTS: list[dict] = [
    {
        "name": "basic_run",
        "script": "examples/basic_run.py",
        "extra_args": ["--plot", str(ASSETS_FIGURES / "basic_run.png")],
    },
    {
        "name": "live_run",
        "script": "examples/live_run.py",
        "extra_args": [],
    },
    {
        "name": "updraft_ensemble",
        "script": "examples/updraft_ensemble.py",
        "extra_args": ["--n", "256", "--plot", str(ASSETS_FIGURES / "updraft_ensemble.png")],
    },
    {
        "name": "differentiable_smax",
        "script": "examples/differentiable_smax.py",
        "extra_args": [],
    },
    {
        "name": "activation_comparison",
        "script": "examples/activation_comparison.py",
        "extra_args": ["--plot", str(ASSETS_FIGURES / "activation_comparison.png")],
    },
    # The sensitivity sweep is expensive (~5 min cold).  It writes a cache to
    # output/sensitivity_sweep_cache.npz and reuses it on subsequent builds if
    # the grid parameters match.
    {
        "name": "sensitivity_sweep",
        "script": "examples/sensitivity_sweep.py",
        "extra_args": ["--plot", str(ASSETS_FIGURES / "sensitivity_sweep.png")],
    },
]


def on_pre_build(config) -> None:  # noqa: ANN001
    ASSETS_FIGURES.mkdir(parents=True, exist_ok=True)
    ASSETS_OUTPUT.mkdir(parents=True, exist_ok=True)

    if not os.environ.get("DOCS_EXECUTE"):
        return

    env = {
        **os.environ,
        "MPLBACKEND": "Agg",
        "JAX_PLATFORM_NAME": "cpu",
    }

    for s in SCRIPTS:
        name = s["name"]
        script = REPO_ROOT / s["script"]
        out_file = ASSETS_OUTPUT / f"{name}.txt"
        cmd = [sys.executable, str(script)] + s["extra_args"]

        print(f"[docs] executing {s['script']} …")
        try:
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                timeout=600,
                env=env,
                cwd=str(REPO_ROOT),
            )
            stdout = result.stdout.strip() or "(no output)"
            out_file.write_text(stdout)
            if result.returncode != 0:
                print(f"[docs] WARNING: {s['script']} exited {result.returncode}")
                if result.stderr:
                    print(result.stderr[-1000:])
        except Exception as exc:
            msg = f"(execution failed: {exc})"
            out_file.write_text(msg)
            print(f"[docs] ERROR running {s['script']}: {exc}")
