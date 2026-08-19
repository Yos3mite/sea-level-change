import numpy as np
import pandas as pd
import pytest
import xarray as xr

from gmsl_budget.gia import apply_piecewise_trend_correction, apply_scalar_gia, area_average_spatial_gia
from gmsl_budget.models import MonthlySeries, SpatialMask


def _series(metadata=None):
    return MonthlySeries(
        pd.to_datetime(["2010-01-15", "2020-01-15"]),
        np.zeros(2),
        "gmsl_raw",
        "mm",
        {"gia_corrected": False, **(metadata or {})},
    )


def test_jin_gia_rate_subtracts_three_mm_over_ten_years():
    raw = _series()
    corrected = apply_scalar_gia(raw, -0.30, reference_time=raw.time[0])
    assert corrected.values[1] - corrected.values[0] == pytest.approx(-3.0, rel=2e-4)
    assert corrected.metadata["gia_corrected"] is True
    assert corrected.metadata["gia_sign_convention"] == "signed rate is added to the parent series"


def test_wet_troposphere_drift_is_piecewise_from_2016():
    raw = MonthlySeries(
        pd.to_datetime(["2015-01-15", "2016-01-15", "2018-01-15"]),
        np.zeros(3),
        "gmsl_raw",
        "mm",
        {},
    )
    corrected = apply_piecewise_trend_correction(
        raw,
        rate_mm_per_year=-0.50,
        start_time=pd.Timestamp("2016-01-15"),
        correction_name="jason3_wet_troposphere_drift",
    )
    assert corrected.values == pytest.approx([0.0, 0.0, -1.0], abs=2e-3)
    assert corrected.metadata["trend_correction_rate_mm_per_year"] == -0.50
    assert corrected.metadata["trend_correction_start_time"] == "2016-01-15T00:00:00"


def test_gia_cannot_be_applied_twice():
    corrected = _series({"gia_corrected": True})
    with pytest.raises(ValueError, match="already GIA corrected"):
        apply_scalar_gia(corrected, 0.30)


def test_official_corrected_indicator_is_rejected_even_if_flag_missing():
    indicator = MonthlySeries(
        pd.to_datetime(["2010-01-15", "2020-01-15"]),
        np.zeros(2),
        "MSL_filtered_GIA_corrected_adjusted",
        "mm",
        {},
    )
    with pytest.raises(ValueError, match="already GIA corrected"):
        apply_scalar_gia(indicator, 0.30)


def test_spatial_gia_average_uses_mask_fraction_and_area():
    rate = xr.DataArray(
        np.array([[0.0, 2.0]]),
        dims=("latitude", "longitude"),
        coords={"latitude": [0.0], "longitude": [0.0, 180.0]},
        attrs={"units": "mm/year"},
    )
    mask = SpatialMask(
        latitude=np.array([0.0]),
        longitude=np.array([0.0, 180.0]),
        ocean_fraction=np.array([[1.0, 0.5]]),
        support=np.array([[True, True]]),
        metadata={},
    )
    assert area_average_spatial_gia(rate, mask) == pytest.approx(2.0 / 3.0)


def test_spatial_gia_rejects_wrong_units():
    rate = xr.DataArray(
        np.zeros((1, 1)),
        dims=("latitude", "longitude"),
        coords={"latitude": [0.0], "longitude": [0.0]},
        attrs={"units": "m/year"},
    )
    mask = SpatialMask(np.array([0.0]), np.array([0.0]), np.ones((1, 1)), np.ones((1, 1), bool), {})
    with pytest.raises(ValueError, match="mm/year"):
        area_average_spatial_gia(rate, mask)
