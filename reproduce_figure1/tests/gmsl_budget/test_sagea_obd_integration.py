from datetime import date
from pathlib import Path

import numpy as np

from gmsl_budget.obd import (
    ProcessedCoefficients,
    convert_coefficients_to_vertical_displacement,
    load_wang_love_numbers,
)


SAGEA_ROOT = Path("D:/AAAA海平面变化/SaGEA-main(L2-L3)/SaGEA-main")


def _processed(c, s):
    return ProcessedCoefficients(
        c=np.asarray([c], dtype=np.float64),
        s=np.asarray([s], dtype=np.float64),
        start=(date(2020, 1, 1),),
        end=(date(2020, 1, 31),),
        midpoint=(date(2020, 1, 16),),
        lmax=c.shape[0] - 1,
        preprocessing_hash="synthetic",
        metadata={},
    )


def test_wang_love_numbers_have_degree_zero_and_requested_length():
    k, h = load_wang_love_numbers(SAGEA_ROOT, lmax=2)
    assert k.shape == (3,)
    assert h.shape == (3,)
    assert k[0] == 0.0
    assert h[0] == 0.0
    assert np.all(np.isfinite(k))
    assert np.all(np.isfinite(h))


def test_zero_coefficients_synthesize_exactly_zero_float64_grid():
    c = np.zeros((3, 3), dtype=np.float64)
    s = np.zeros_like(c)
    grid = convert_coefficients_to_vertical_displacement(
        _processed(c, s),
        SAGEA_ROOT,
        latitude=np.array([-45.0, 45.0]),
        longitude=np.array([0.0, 180.0]),
    )
    assert grid.dtype == np.float64
    assert grid.shape == (1, 2, 2)
    assert np.array_equal(grid.values, np.zeros((1, 2, 2)))
    assert grid.attrs["preprocessing_hash"] == "synthetic"
