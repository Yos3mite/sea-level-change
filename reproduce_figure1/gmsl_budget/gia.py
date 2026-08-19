from __future__ import annotations

import numpy as np
import pandas as pd
import xarray as xr

from .masks import cell_area_weights
from .models import MonthlySeries, SpatialMask
from .trend import decimal_year


def _is_already_corrected(series: MonthlySeries) -> bool:
    if series.metadata.get("gia_corrected") is True:
        return True
    labels = [series.name, str(series.metadata.get("source_variable", ""))]
    return any("gia_corrected" in label.lower() for label in labels)


def apply_scalar_gia(
    series: MonthlySeries,
    rate_mm_per_year: float = 0.30,
    reference_time: pd.Timestamp | None = None,
) -> MonthlySeries:
    if series.units != "mm":
        raise ValueError("scalar GIA correction requires millimetres (mm)")
    if _is_already_corrected(series):
        raise ValueError("input series is already GIA corrected")
    if not np.isfinite(rate_mm_per_year):
        raise ValueError("GIA rate must be finite")
    reference = pd.Timestamp(series.time[0] if reference_time is None else reference_time)
    offset_years = decimal_year(series.time) - decimal_year(pd.DatetimeIndex([reference]))[0]
    corrected = series.values + float(rate_mm_per_year) * offset_years
    metadata = dict(series.metadata)
    metadata.update(
        {
            "gia_corrected": True,
            "gia_mode": "scalar",
            "gia_rate_mm_per_year": float(rate_mm_per_year),
            "gia_reference_time": reference.isoformat(),
            "gia_sign_convention": "positive rate increases corrected GMSL trend",
            "parent_series": series.name,
        }
    )
    return MonthlySeries(series.time, corrected, "gmsl_gia_corrected", "mm", metadata)


def area_average_spatial_gia(gia_rate_grid: xr.DataArray, mask: SpatialMask) -> float:
    units = str(gia_rate_grid.attrs.get("units", "")).strip().lower()
    if units not in {"mm/year", "mm/yr", "mm yr-1", "mm a-1"}:
        raise ValueError("spatial GIA rate units must be mm/year")
    grid = gia_rate_grid.transpose("latitude", "longitude")
    if not np.allclose(grid.latitude.values, mask.latitude):
        raise ValueError("spatial GIA latitude does not match mask")
    if not np.allclose(grid.longitude.values, mask.longitude):
        raise ValueError("spatial GIA longitude does not match mask")
    values = np.asarray(grid.values, dtype=np.float64)
    weights = cell_area_weights(mask.latitude, mask.longitude) * mask.ocean_fraction * mask.support
    valid = np.isfinite(values) & mask.support
    valid_weight = float(np.sum(weights[valid]))
    total_weight = float(np.sum(weights))
    if total_weight <= 0 or valid_weight / total_weight < 0.995:
        raise ValueError("spatial GIA grid does not cover at least 99.5% of the fixed mask")
    return float(np.sum(np.where(valid, values, 0.0) * weights) / total_weight)
