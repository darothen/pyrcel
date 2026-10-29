"""Standalone regressions for this fix on unmodified upstream."""

import numpy as np
import pytest

from pyrcel import Lognorm


def test_distribution_moment_identities():
    d = Lognorm(0.05, 2.0, 100.0)
    s = d.stats()
    assert s["effective_radius"] == pytest.approx(0.05 * np.exp(2.5 * np.log(2.0) ** 2))
    assert s["total_diameter"] == pytest.approx(2 * 100 * 0.05 * np.exp(0.5 * np.log(2.0) ** 2))
