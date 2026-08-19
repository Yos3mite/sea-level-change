from datetime import date

import numpy as np
import pandas as pd
import pytest
import xarray as xr

from gmsl_budget.models import SpatialMask
from gmsl_budget.obd import (
    EARTH_RADIUS_M,
    ProcessedCoefficients,
    area_average_obd,
    assert_shared_preprocessing,
    vertical_conversion_weights,
)


def test_vertical_displacement_uses_h_over_one_plus_k():
    k = np.array([0.0, 0.0, -0.3])
    h = np.array([0.0, 0.0, -0.6])
    got = vertical_conversion_weights(lmax=2, k=k, h=h)
    assert got[2] == pytest.approx(EARTH_RADIUS_M * (-0.6) / 0.7)


def test_vertical_weights_reject_singular_one_plus_k():
    with pytest.raises(ValueError, match=r"1 \+ k"):
        vertical_conversion_weights(1, np.array([0.0, -1.0]), np.array([0.0, 0.1]))


def test_downward_ocean_bottom_motion_is_negative_obd():
    displacement = xr.DataArray(
        np.array([[[-0.002, -0.002]]]),
        dims=("time", "latitude", "longitude"),
        coords={"time": pd.to_datetime(["2020-01-15"]), "latitude": [0.0], "longitude": [0.0, 180.0]},
        attrs={"units": "m", "preprocessing_hash": "abc123"},
    )
    mask = SpatialMask(
        np.array([0.0]),
        np.array([0.0, 180.0]),
        np.ones((1, 2)),
        np.ones((1, 2), dtype=bool),
        {"sha256": "mask"},
    )
    result = area_average_obd(displacement, mask, "abc123")
    assert result.values.tolist() == pytest.approx([-2.0])
    assert result.metadata["sign_convention"] == "upward positive; subsidence negative"


def test_mass_and_obd_reject_different_preprocessing_hashes():
    with pytest.raises(ValueError, match="preprocessing hash"):
        assert_shared_preprocessing(mass_hash="aaa", obd_hash="bbb")


def test_processed_coefficients_reject_duplicate_months():
    c = np.zeros((2, 3, 3))
    s = np.zeros_like(c)
    with pytest.raises(ValueError, match="duplicate month"):
        ProcessedCoefficients(
            c=c,
            s=s,
            start=(date(2020, 1, 1), date(2020, 1, 10)),
            end=(date(2020, 1, 31), date(2020, 1, 31)),
            midpoint=(date(2020, 1, 16), date(2020, 1, 20)),
            lmax=2,
            preprocessing_hash="hash",
            metadata={},
        )
