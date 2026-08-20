"""Validated common data types for regional Figure 3 processing."""

from __future__ import annotations

from dataclasses import dataclass
import re
from typing import Any, Mapping

import numpy as np


_MONTH_PATTERN = re.compile(r"^\d{4}-(0[1-9]|1[0-2])$")
_MONTH_STATUSES = frozenset({"observed", "reconstructed", "missing"})


@dataclass(frozen=True)
class MonthlyGridSeries:
    """One monthly EWH grid product normalized to ``(time, lat, lon)``."""

    source_id: str
    months: np.ndarray
    lat: np.ndarray
    lon: np.ndarray
    ewh_mm: np.ndarray
    valid_month: np.ndarray
    month_status: np.ndarray
    metadata: Mapping[str, Any]

    def __post_init__(self) -> None:
        months = np.asarray(self.months, dtype="U7")
        lat = np.asarray(self.lat, dtype=float)
        lon = np.asarray(self.lon, dtype=float)
        ewh_mm = np.asarray(self.ewh_mm, dtype=float)
        valid_month = np.asarray(self.valid_month, dtype=bool)
        month_status = np.asarray(self.month_status, dtype="U13")

        if not self.source_id:
            raise ValueError("source_id must be non-empty")
        if months.ndim != 1 or any(not _MONTH_PATTERN.fullmatch(x) for x in months):
            raise ValueError("months must be a one-dimensional YYYY-MM array")
        if len(np.unique(months)) != len(months):
            raise ValueError("duplicate months are not allowed")
        if months.size > 1 and np.any(months[1:] <= months[:-1]):
            raise ValueError("months must be strictly increasing")

        for name, coordinate in (("lat", lat), ("lon", lon)):
            if coordinate.ndim != 1 or not np.isfinite(coordinate).all():
                raise ValueError(f"{name} must be a finite one-dimensional array")
            if coordinate.size > 1 and np.any(np.diff(coordinate) <= 0):
                raise ValueError(f"{name} must be strictly increasing")

        expected_shape = (months.size, lat.size, lon.size)
        if ewh_mm.shape != expected_shape:
            raise ValueError(
                f"ewh_mm shape must be {expected_shape}, got {ewh_mm.shape}"
            )
        if valid_month.shape != (months.size,):
            raise ValueError("valid_month must contain one value per month")
        if month_status.shape != (months.size,):
            raise ValueError("month_status must contain one value per month")
        unknown = sorted(set(month_status.tolist()) - _MONTH_STATUSES)
        if unknown:
            raise ValueError(f"unsupported month status: {unknown}")
        if np.any((month_status == "missing") != ~valid_month):
            raise ValueError("missing status and valid_month flags are inconsistent")
        if np.any([not np.isnan(ewh_mm[i]).all() for i in np.flatnonzero(~valid_month)]):
            raise ValueError("missing months must contain only NaN fields")

        object.__setattr__(self, "months", months)
        object.__setattr__(self, "lat", lat)
        object.__setattr__(self, "lon", lon)
        object.__setattr__(self, "ewh_mm", ewh_mm)
        object.__setattr__(self, "valid_month", valid_month)
        object.__setattr__(self, "month_status", month_status)
        object.__setattr__(self, "metadata", dict(self.metadata))

