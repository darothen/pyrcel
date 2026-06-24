"""Smoke tests for the run_parcel CLI (issue #60)."""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path

import pytest

FIXTURE = Path(__file__).parent / "fixtures" / "cli_simple.yml"


@pytest.mark.slow
def test_cli_runs_and_writes_netcdf(tmp_path):
    """run_parcel with a minimal YAML produces a valid NetCDF file."""
    out = tmp_path / "smoke.nc"
    result = subprocess.run(
        [sys.executable, "-m", "pyrcel.scripts.run_parcel", str(FIXTURE), "-o", str(out),
         "--no-console"],
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0, f"CLI failed:\n{result.stderr}"
    assert out.exists(), "NetCDF file was not created"
    assert out.stat().st_size > 0

    import xarray as xr

    ds = xr.open_dataset(out)
    assert "S" in ds
    assert "T" in ds
    assert ds.sizes["time"] > 1
    ds.close()


@pytest.mark.slow
def test_cli_output_dir_from_namelist(tmp_path):
    """run_parcel writes to <output_dir>/<name>.nc when -o is omitted."""
    import tempfile

    import yaml

    cfg = yaml.safe_load(FIXTURE.read_text())
    cfg["experiment_control"]["output_dir"] = str(tmp_path / "out")
    cfg["experiment_control"]["name"] = "test_run"

    with tempfile.NamedTemporaryFile(suffix=".yml", mode="w", delete=False) as f:
        yaml.dump(cfg, f)
        yml_path = f.name

    try:
        result = subprocess.run(
            [sys.executable, "-m", "pyrcel.scripts.run_parcel", yml_path, "--no-console"],
            capture_output=True,
            text=True,
        )
        assert result.returncode == 0, f"CLI failed:\n{result.stderr}"
        expected = tmp_path / "out" / "test_run.nc"
        assert expected.exists()
    finally:
        Path(yml_path).unlink(missing_ok=True)


def test_cli_missing_namelist(tmp_path):
    """run_parcel exits 1 when the namelist file does not exist."""
    result = subprocess.run(
        [sys.executable, "-m", "pyrcel.scripts.run_parcel", str(tmp_path / "nope.yml")],
        capture_output=True,
        text=True,
    )
    assert result.returncode == 1


def test_cli_help():
    """run_parcel --help exits 0 and mentions the namelist argument."""
    result = subprocess.run(
        [sys.executable, "-m", "pyrcel.scripts.run_parcel", "--help"],
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0
    assert "config.yml" in result.stdout
