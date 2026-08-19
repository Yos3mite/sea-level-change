import numpy as np
import pandas as pd
import pytest
from scipy.io import savemat
import xarray as xr

from gmsl_budget.cmems import build_altimetry_mask, read_cmems_gmsl
from gmsl_budget.masks import (
    EARTH_RADIUS_M,
    buffer_ocean_mask,
    cell_area_weights,
    load_cdt_coast_distance,
    load_cdt_ocean_fraction,
)
from gmsl_budget.models import SpatialMask


def test_area_weights_use_cell_bounds_not_only_cosine_centers():
    latitude = np.array([-60.0, 0.0, 60.0])
    longitude = np.array([0.0, 120.0, 240.0])
    weights = cell_area_weights(latitude, longitude)
    assert weights[1, 0] > weights[0, 0]
    assert weights[:, 0].sum() == pytest.approx(4 * np.pi * EARTH_RADIUS_M**2 / 3, rel=1e-12)


def test_monthly_missing_cell_does_not_renormalize_fixed_support(tmp_path):
    path = tmp_path / "cmems.nc"
    data = xr.Dataset(
        {
            "sla": (
                ("time", "latitude", "longitude"),
                np.array([[[1.0, 3.0]], [[1.0, np.nan]]], dtype=np.float64),
                {"units": "m"},
            )
        },
        coords={
            "time": pd.to_datetime(["2020-01-01", "2020-02-01"]),
            "latitude": np.array([0.0]),
            "longitude": np.array([0.0, 180.0]),
        },
    )
    data.to_netcdf(path, engine="scipy")
    mask = SpatialMask(
        latitude=np.array([0.0]),
        longitude=np.array([0.0, 180.0]),
        ocean_fraction=np.ones((1, 2)),
        support=np.ones((1, 2), dtype=bool),
        metadata={"name": "fixed"},
    )
    result = read_cmems_gmsl(path, mask, min_valid_weight_fraction=0.995)
    assert result.values[0] == pytest.approx(2000.0)
    assert np.isnan(result.values[1])
    assert result.metadata["valid_weight_fraction"] == pytest.approx([1.0, 0.5])


def test_build_altimetry_mask_uses_common_valid_cells_across_time():
    sla = xr.DataArray(
        np.array([[[1.0, 2.0]], [[1.0, np.nan]]]),
        dims=("time", "latitude", "longitude"),
        coords={"time": [0, 1], "latitude": [0.0], "longitude": [0.0, 180.0]},
    )
    fraction = xr.DataArray(
        np.array([[1.0, 0.5]]),
        dims=("latitude", "longitude"),
        coords={"latitude": [0.0], "longitude": [0.0, 180.0]},
    )
    mask = build_altimetry_mask(sla, fraction)
    assert mask.support.tolist() == [[True, False]]
    assert mask.ocean_fraction.tolist() == [[1.0, 0.5]]


def test_buffer_keeps_only_ocean_cells_at_least_300_km_from_coast():
    mask = SpatialMask(
        latitude=np.array([0.0]),
        longitude=np.array([0.0, 1.0, 2.0]),
        ocean_fraction=np.ones((1, 3)),
        support=np.ones((1, 3), dtype=bool),
        metadata={},
    )
    distance = xr.DataArray(
        np.array([[150.0, 300.0, 350.0]]),
        dims=("latitude", "longitude"),
        coords={"latitude": [0.0], "longitude": [0.0, 1.0, 2.0]},
    )
    buffered = buffer_ocean_mask(mask, distance, distance_km=300.0)
    assert buffered.support.tolist() == [[False, True, True]]
    assert buffered.metadata["coastal_buffer_km"] == 300.0


def test_cdt_subcell_aggregation_returns_fractional_coastal_cell(tmp_path):
    land_path = tmp_path / "land_mask.mat"
    savemat(
        land_path,
        {
            "lat": np.array([[-0.0625], [0.0625]]),
            "lon": np.array([[-0.0625, 0.0625]]),
            "land": np.array([[0, 1], [0, 1]], dtype=np.uint8),
        },
    )
    fraction = load_cdt_ocean_fraction(
        land_path,
        target_latitude=np.array([0.0]),
        target_longitude=np.array([0.0]),
        target_spacing_degrees=0.25,
    )
    assert fraction.item() == pytest.approx(0.5)


def test_cdt_distance_loader_samples_target_centres(tmp_path):
    distance_path = tmp_path / "distance2coast.mat"
    savemat(
        distance_path,
        {
            "lat": np.array([[-1.0], [0.0], [1.0]]),
            "lon": np.array([[-1.0, 0.0, 1.0]]),
            "D": np.array([[100.0, 200.0, 300.0], [400.0, 500.0, 600.0], [700.0, 800.0, 900.0]]),
        },
    )
    result = load_cdt_coast_distance(
        distance_path,
        target_latitude=np.array([0.0]),
        target_longitude=np.array([0.0]),
    )
    assert result.item() == 500.0
