from __future__ import annotations

from pathlib import Path

import numpy as np
from scipy.io import loadmat
import xarray as xr

from .models import SpatialMask


EARTH_RADIUS_M = 6_371_000.0
CDT_CITATION = "Greene et al. (2019), doi:10.1029/2019GC008392"


def _latitude_bounds(latitude: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    latitude = np.asarray(latitude, dtype=np.float64)
    if latitude.ndim != 1 or len(latitude) == 0:
        raise ValueError("latitude must be a non-empty one-dimensional array")
    order = np.argsort(latitude)
    sorted_latitude = latitude[order]
    if len(np.unique(sorted_latitude)) != len(sorted_latitude):
        raise ValueError("latitude must be unique")
    if np.any((sorted_latitude < -90.0) | (sorted_latitude > 90.0)):
        raise ValueError("latitude must be within [-90, 90]")
    if len(sorted_latitude) == 1:
        lower_sorted = np.array([-90.0])
        upper_sorted = np.array([90.0])
    else:
        middle = (sorted_latitude[:-1] + sorted_latitude[1:]) / 2.0
        lower_sorted = np.concatenate(([-90.0], middle))
        upper_sorted = np.concatenate((middle, [90.0]))
    lower = np.empty_like(lower_sorted)
    upper = np.empty_like(upper_sorted)
    lower[order] = lower_sorted
    upper[order] = upper_sorted
    return lower, upper


def _longitude_widths(longitude: np.ndarray) -> np.ndarray:
    longitude = np.asarray(longitude, dtype=np.float64)
    if longitude.ndim != 1 or len(longitude) == 0:
        raise ValueError("longitude must be a non-empty one-dimensional array")
    normalized = np.mod(longitude, 360.0)
    if len(np.unique(np.round(normalized, 12))) != len(normalized):
        raise ValueError("longitude must be unique modulo 360 degrees")
    if len(longitude) == 1:
        return np.array([360.0])
    order = np.argsort(normalized)
    sorted_longitude = normalized[order]
    forward_gaps = np.diff(np.concatenate((sorted_longitude, [sorted_longitude[0] + 360.0])))
    widths_sorted = (np.roll(forward_gaps, 1) + forward_gaps) / 2.0
    widths = np.empty_like(widths_sorted)
    widths[order] = widths_sorted
    return widths


def cell_area_weights(latitude: np.ndarray, longitude: np.ndarray) -> np.ndarray:
    lower, upper = _latitude_bounds(latitude)
    longitude_width = np.deg2rad(_longitude_widths(longitude))
    latitude_factor = np.sin(np.deg2rad(upper)) - np.sin(np.deg2rad(lower))
    return EARTH_RADIUS_M**2 * latitude_factor[:, None] * longitude_width[None, :]


def _nearest_latitude_indices(source: np.ndarray, query: np.ndarray) -> np.ndarray:
    source = np.asarray(source, dtype=np.float64).reshape(-1)
    return np.asarray([int(np.argmin(np.abs(source - value))) for value in query], dtype=int)


def _nearest_longitude_indices(source: np.ndarray, query: np.ndarray) -> np.ndarray:
    source = np.asarray(source, dtype=np.float64).reshape(-1)
    indices = []
    for value in query:
        distance = np.abs((source - value + 180.0) % 360.0 - 180.0)
        indices.append(int(np.argmin(distance)))
    return np.asarray(indices, dtype=int)


def load_cdt_ocean_fraction(
    land_mask_mat: str | Path,
    target_latitude: np.ndarray,
    target_longitude: np.ndarray,
    target_spacing_degrees: float = 0.25,
) -> xr.DataArray:
    if target_spacing_degrees <= 0:
        raise ValueError("target_spacing_degrees must be positive")
    source = loadmat(Path(land_mask_mat))
    source_latitude = np.asarray(source["lat"], dtype=np.float64).reshape(-1)
    source_longitude = np.asarray(source["lon"], dtype=np.float64).reshape(-1)
    land = np.asarray(source["land"], dtype=np.float64)
    if land.shape != (len(source_latitude), len(source_longitude)):
        raise ValueError("CDT land mask coordinates do not match land array")
    target_latitude = np.asarray(target_latitude, dtype=np.float64)
    target_longitude = np.asarray(target_longitude, dtype=np.float64)
    offset = target_spacing_degrees / 4.0
    latitude_samples = np.column_stack((target_latitude - offset, target_latitude + offset)).reshape(-1)
    longitude_samples = np.column_stack((target_longitude - offset, target_longitude + offset)).reshape(-1)
    lat_index = _nearest_latitude_indices(source_latitude, latitude_samples)
    lon_index = _nearest_longitude_indices(source_longitude, longitude_samples)
    sampled_land = land[np.ix_(lat_index, lon_index)].reshape(
        len(target_latitude), 2, len(target_longitude), 2
    )
    ocean_fraction = 1.0 - sampled_land.mean(axis=(1, 3))
    return xr.DataArray(
        ocean_fraction.astype(np.float64),
        dims=("latitude", "longitude"),
        coords={"latitude": target_latitude, "longitude": target_longitude},
        attrs={"source": str(Path(land_mask_mat).resolve()), "citation": CDT_CITATION},
    )


def load_cdt_coast_distance(
    distance_mat: str | Path,
    target_latitude: np.ndarray,
    target_longitude: np.ndarray,
) -> xr.DataArray:
    source = loadmat(Path(distance_mat))
    source_latitude = np.asarray(source["lat"], dtype=np.float64).reshape(-1)
    source_longitude = np.asarray(source["lon"], dtype=np.float64).reshape(-1)
    distance = np.asarray(source["D"], dtype=np.float64)
    if distance.shape != (len(source_latitude), len(source_longitude)):
        raise ValueError("CDT distance coordinates do not match D array")
    target_latitude = np.asarray(target_latitude, dtype=np.float64)
    target_longitude = np.asarray(target_longitude, dtype=np.float64)
    lat_index = _nearest_latitude_indices(source_latitude, target_latitude)
    lon_index = _nearest_longitude_indices(source_longitude, target_longitude)
    sampled = distance[np.ix_(lat_index, lon_index)]
    return xr.DataArray(
        sampled.astype(np.float64),
        dims=("latitude", "longitude"),
        coords={"latitude": target_latitude, "longitude": target_longitude},
        attrs={
            "units": "km",
            "distance_method": "Haversine great-circle distance",
            "source": str(Path(distance_mat).resolve()),
            "citation": CDT_CITATION,
        },
    )


def buffer_ocean_mask(
    mask: SpatialMask,
    coast_distance_km: xr.DataArray,
    distance_km: float,
) -> SpatialMask:
    if distance_km < 0:
        raise ValueError("distance_km must be non-negative")
    distance = np.asarray(coast_distance_km.values, dtype=np.float64)
    if distance.shape != mask.support.shape:
        raise ValueError("coast distance grid does not match mask shape")
    if not np.allclose(coast_distance_km.latitude.values, mask.latitude):
        raise ValueError("coast distance latitude does not match mask")
    if not np.allclose(coast_distance_km.longitude.values, mask.longitude):
        raise ValueError("coast distance longitude does not match mask")
    support = mask.support & (mask.ocean_fraction > 0.0) & np.isfinite(distance) & (distance >= distance_km)
    metadata = dict(mask.metadata)
    metadata.update(
        {
            "coastal_buffer_km": float(distance_km),
            "coast_distance_source": coast_distance_km.attrs.get("source"),
            "coast_distance_method": coast_distance_km.attrs.get("distance_method"),
        }
    )
    return SpatialMask(mask.latitude, mask.longitude, mask.ocean_fraction, support, metadata)
