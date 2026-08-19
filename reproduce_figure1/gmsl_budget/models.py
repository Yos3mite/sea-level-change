from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Mapping

import numpy as np
import pandas as pd


def _readonly_float64(values: np.ndarray) -> np.ndarray:
    array = np.array(values, dtype=np.float64, copy=True)
    array.setflags(write=False)
    return array


@dataclass(frozen=True)
class MonthlySeries:
    time: pd.DatetimeIndex
    values: np.ndarray
    name: str
    units: str
    metadata: Mapping[str, Any]

    def __post_init__(self) -> None:
        time = pd.DatetimeIndex(pd.to_datetime(self.time))
        if time.hasnans:
            raise ValueError("time contains NaT")
        months = time.to_period("M")
        if months.duplicated().any():
            raise ValueError("duplicate month in MonthlySeries")
        normalized_time = pd.DatetimeIndex(
            months.to_timestamp(how="start") + pd.Timedelta(days=14),
        ).rename(None)
        values = _readonly_float64(self.values)
        if values.shape != (len(normalized_time),):
            raise ValueError("values must have shape (time,)")
        if not self.name:
            raise ValueError("name must not be empty")
        if not self.units:
            raise ValueError("units must not be empty")
        object.__setattr__(self, "time", normalized_time)
        object.__setattr__(self, "values", values)
        object.__setattr__(self, "metadata", dict(self.metadata))


@dataclass(frozen=True)
class SpatialMask:
    latitude: np.ndarray
    longitude: np.ndarray
    ocean_fraction: np.ndarray
    support: np.ndarray
    metadata: Mapping[str, Any]

    def __post_init__(self) -> None:
        latitude = _readonly_float64(self.latitude)
        longitude = _readonly_float64(self.longitude)
        fraction = _readonly_float64(self.ocean_fraction)
        support = np.array(self.support, dtype=bool, copy=True)
        support.setflags(write=False)
        if latitude.ndim != 1 or longitude.ndim != 1:
            raise ValueError("latitude and longitude must be one-dimensional")
        expected_shape = (len(latitude), len(longitude))
        if fraction.shape != expected_shape or support.shape != expected_shape:
            raise ValueError(f"mask arrays must have shape {expected_shape}")
        if not np.all(np.isfinite(latitude)) or not np.all(np.isfinite(longitude)):
            raise ValueError("coordinates must be finite")
        if len(np.unique(latitude)) != len(latitude) or len(np.unique(longitude)) != len(longitude):
            raise ValueError("coordinates must be unique")
        if np.any(~np.isfinite(fraction)) or np.any((fraction < 0.0) | (fraction > 1.0)):
            raise ValueError("ocean_fraction must be finite and within [0, 1]")
        object.__setattr__(self, "latitude", latitude)
        object.__setattr__(self, "longitude", longitude)
        object.__setattr__(self, "ocean_fraction", fraction)
        object.__setattr__(self, "support", support)
        object.__setattr__(self, "metadata", dict(self.metadata))


@dataclass(frozen=True)
class TrendResult:
    series_name: str
    trend_mm_per_year: float
    ols_standard_error: float
    hac_standard_error: float
    intercept_mm: float
    annual_sin_mm: float
    annual_cos_mm: float
    semiannual_sin_mm: float
    semiannual_cos_mm: float
    residual_lag1_correlation: float
    n_obs: int
    start_time: pd.Timestamp
    end_time: pd.Timestamp
    missing_months: tuple[str, ...]
