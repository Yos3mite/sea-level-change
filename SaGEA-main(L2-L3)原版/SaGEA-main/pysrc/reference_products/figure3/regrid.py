"""Deterministic regridding for normalized Figure 3 monthly products."""

from __future__ import annotations

import numpy as np
from scipy.interpolate import RegularGridInterpolator

from .types import MonthlyGridSeries


def nearest_regrid(
    series: MonthlyGridSeries,
    target_lat: np.ndarray,
    target_lon: np.ndarray,
) -> MonthlyGridSeries:
    """Nearest-neighbor regrid with a cyclic longitude seam."""
    target_lat = np.asarray(target_lat, dtype=float)
    target_lon = np.asarray(target_lon, dtype=float)
    if target_lat.ndim != 1 or target_lon.ndim != 1:
        raise ValueError("target latitude and longitude must be one-dimensional")
    if np.any(np.diff(target_lat) <= 0) or np.any(np.diff(target_lon) <= 0):
        raise ValueError("target latitude and longitude must be strictly increasing")

    normalized_target_lon = (target_lon + 180.0) % 360.0 - 180.0
    if np.any(np.diff(normalized_target_lon) <= 0):
        raise ValueError("target longitude must remain ordered after normalization")
    cyclic_lon = np.concatenate(
        ([series.lon[-1] - 360.0], series.lon, [series.lon[0] + 360.0])
    )
    lat_mesh, lon_mesh = np.meshgrid(
        target_lat, normalized_target_lon, indexing="ij"
    )
    query = np.column_stack((lat_mesh.ravel(), lon_mesh.ravel()))
    output = np.full(
        (series.months.size, target_lat.size, target_lon.size), np.nan, dtype=float
    )

    for month_index in range(series.months.size):
        field = series.ewh_mm[month_index]
        cyclic_field = np.concatenate(
            (field[:, -1:], field, field[:, :1]), axis=1
        )
        interpolator = RegularGridInterpolator(
            (series.lat, cyclic_lon),
            cyclic_field,
            method="nearest",
            bounds_error=False,
            fill_value=np.nan,
        )
        output[month_index] = interpolator(query).reshape(
            target_lat.size, target_lon.size
        )

    metadata = dict(series.metadata)
    metadata["regridding"] = {
        "method": "nearest",
        "longitude_cyclic": True,
        "target_shape": [int(target_lat.size), int(target_lon.size)],
    }
    return MonthlyGridSeries(
        source_id=series.source_id,
        months=series.months,
        lat=target_lat,
        lon=normalized_target_lon,
        ewh_mm=output,
        valid_month=series.valid_month,
        month_status=series.month_status,
        metadata=metadata,
    )
