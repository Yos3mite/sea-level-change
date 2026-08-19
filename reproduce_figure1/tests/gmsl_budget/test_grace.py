import numpy as np
import pandas as pd
import pytest
import xarray as xr

from gmsl_budget.grace import (
    align_common_months,
    ensemble_mean,
    ewh_to_sea_level,
    read_mascon_ocean_series,
)
from gmsl_budget.models import MonthlySeries, SpatialMask


def _series(months, values, center):
    return MonthlySeries(
        pd.PeriodIndex(months, freq="M").to_timestamp() + pd.Timedelta(days=14),
        np.asarray(values, dtype=np.float64),
        f"ocean_mass_{center.lower()}",
        "mm",
        {"center": center, "preprocessing_hash": "same", "mask_hash": "mask"},
    )


def test_ewh_density_conversion_uses_fresh_to_seawater_ratio():
    assert ewh_to_sea_level(np.array([1028.0]))[0] == pytest.approx(1000.0)


def test_alignment_uses_year_month_keys_not_array_positions():
    first = _series(["2020-01", "2020-03"], [1.0, 3.0], "CSR")
    second = _series(["2020-02", "2020-03"], [2.0, 4.0], "JPL")
    aligned_first, aligned_second = align_common_months([first, second])
    assert list(aligned_first.time.to_period("M").astype(str)) == ["2020-03"]
    assert aligned_first.values.tolist() == [3.0]
    assert aligned_second.values.tolist() == [4.0]


def test_ensemble_metadata_preserves_actual_centers():
    result = ensemble_mean(
        {
            "JPL": _series(["2020-01", "2020-02"], [2.0, 4.0], "JPL"),
            "CSR": _series(["2020-01", "2020-02"], [0.0, 2.0], "CSR"),
            "GFZ": _series(["2020-01", "2020-02"], [1.0, 3.0], "GFZ"),
        }
    )
    assert result.values.tolist() == pytest.approx([1.0, 3.0])
    assert result.metadata["centers"] == ["CSR", "GFZ", "JPL"]
    assert result.metadata["ensemble_spread_mm"] == pytest.approx([1.0, 1.0])


def test_ensemble_rejects_different_spatial_masks():
    csr = _series(["2020-01"], [0.0], "CSR")
    jpl = MonthlySeries(
        csr.time,
        np.array([1.0]),
        "ocean_mass_jpl",
        "mm",
        {"center": "JPL", "preprocessing_hash": "different", "mask_hash": "other-mask"},
    )
    with pytest.raises(ValueError, match="mask_hash"):
        ensemble_mean({"CSR": csr, "JPL": jpl})


def test_read_mascon_rejects_center_label_disagreement(tmp_path):
    path = tmp_path / "jpl.nc"
    xr.Dataset(
        {
            "lwe_thickness": (
                ("time", "latitude", "longitude"),
                np.ones((1, 1, 1)),
                {"units": "cm"},
            )
        },
        coords={"time": pd.to_datetime(["2020-01-15"]), "latitude": [0.0], "longitude": [0.0]},
        attrs={"institution": "NASA/JPL", "source": "JPL RL06.3"},
    ).to_netcdf(path, engine="scipy")
    mask = SpatialMask(np.array([0.0]), np.array([0.0]), np.ones((1, 1)), np.ones((1, 1), bool), {})
    with pytest.raises(ValueError, match="does not match"):
        read_mascon_ocean_series(path, "CSR", mask)


def test_read_mascon_converts_cm_ewh_to_mm_sea_level(tmp_path):
    path = tmp_path / "csr.nc"
    xr.Dataset(
        {
            "lwe_thickness": (
                ("time", "latitude", "longitude"),
                np.array([[[1.028]], [[2.056]]]),
                {"units": "cm"},
            )
        },
        coords={
            "time": pd.to_datetime(["2020-01-15", "2020-02-15"]),
            "latitude": [0.0],
            "longitude": [0.0],
        },
        attrs={"institution": "Center for Space Research (CSR)", "source": "CSR RL06.3"},
    ).to_netcdf(path, engine="scipy")
    mask = SpatialMask(np.array([0.0]), np.array([0.0]), np.ones((1, 1)), np.ones((1, 1), bool), {})
    result = read_mascon_ocean_series(path, "CSR", mask)
    assert result.values.tolist() == pytest.approx([10.0, 20.0])
    assert result.metadata["rho_freshwater_kg_m3"] == 1000.0
    assert result.metadata["rho_seawater_kg_m3"] == 1028.0
