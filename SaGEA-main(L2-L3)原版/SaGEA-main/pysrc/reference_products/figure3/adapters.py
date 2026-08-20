"""Configuration-driven adapters for Mascon and custom Level-3 grids."""

from __future__ import annotations

from pathlib import Path
from typing import Any, Mapping

from netCDF4 import Dataset, num2date
import numpy as np

from .types import MonthlyGridSeries


_UNIT_SCALE_TO_MM = {
    "m": 1000.0,
    "meter": 1000.0,
    "meters": 1000.0,
    "metre": 1000.0,
    "metres": 1000.0,
    "cm": 10.0,
    "centimeter": 10.0,
    "centimeters": 10.0,
    "centimetre": 10.0,
    "centimetres": 10.0,
    "mm": 1.0,
    "millimeter": 1.0,
    "millimeters": 1.0,
    "millimetre": 1.0,
    "millimetres": 1.0,
}


def _attribute(variable: Any, *names: str, default: Any = None) -> Any:
    for name in names:
        if name in variable.ncattrs():
            return variable.getncattr(name)
    return default


def _as_float_array(variable: Any) -> np.ndarray:
    return np.asarray(np.ma.filled(variable[:], np.nan), dtype=float)


def _load_grid(config: Mapping[str, Any], *, default_status: str) -> MonthlyGridSeries:
    path = Path(config["path"])
    variables = config["variables"]
    dimensions = config["dimensions"]

    with Dataset(path) as dataset:
        time_variable = dataset.variables[variables["time"]]
        time_units = _attribute(time_variable, "units", "Units")
        if not time_units:
            raise ValueError("time variable is missing units")
        calendar = _attribute(time_variable, "calendar", "Calendar", default="standard")
        decoded = num2date(time_variable[:], units=time_units, calendar=calendar)
        months = np.asarray([f"{item.year:04d}-{item.month:02d}" for item in decoded])

        lat = _as_float_array(dataset.variables[variables["lat"]])
        lon = _as_float_array(dataset.variables[variables["lon"]])
        if lat.ndim != 1 or lon.ndim != 1:
            raise ValueError("only one-dimensional latitude/longitude grids are supported")

        field_variable = dataset.variables[variables["field"]]
        field = _as_float_array(field_variable)
        desired_dimensions = (
            dimensions["time"],
            dimensions["lat"],
            dimensions["lon"],
        )
        try:
            transpose_order = tuple(
                field_variable.dimensions.index(name) for name in desired_dimensions
            )
        except ValueError as exc:
            raise ValueError(
                f"field dimensions {field_variable.dimensions} do not contain "
                f"configured dimensions {desired_dimensions}"
            ) from exc
        if field.ndim != 3:
            raise ValueError("field variable must have exactly three dimensions")
        field = np.transpose(field, transpose_order)

        source_units = str(
            config.get("units")
            or _attribute(field_variable, "units", "Units", default="")
        ).strip()
        unit_key = source_units.lower()
        if unit_key not in _UNIT_SCALE_TO_MM:
            raise ValueError(f"unsupported water-height units: {source_units!r}")
        field *= _UNIT_SCALE_TO_MM[unit_key]

        lat_order = np.argsort(lat)
        normalized_lon = (lon + 180.0) % 360.0 - 180.0
        lon_order = np.argsort(normalized_lon)
        lat = lat[lat_order]
        lon = normalized_lon[lon_order]
        field = field[:, lat_order, :][:, :, lon_order]

        valid_name = variables.get("valid_month")
        if valid_name:
            valid_values = np.ma.filled(dataset.variables[valid_name][:], 0)
            valid_month = np.asarray(valid_values, dtype=bool)
        else:
            valid_month = np.any(np.isfinite(field), axis=(1, 2))

    if valid_month.shape != (months.size,):
        raise ValueError("valid_month variable must have the time dimension only")
    field[~valid_month, :, :] = np.nan
    status = str(config.get("month_status", default_status))
    month_status = np.full(months.size, status, dtype="U13")
    month_status[~valid_month] = "missing"

    metadata = dict(config.get("metadata", {}))
    metadata.update(
        {
            "source_path": str(path),
            "source_units": source_units,
            "normalized_units": "mm",
        }
    )
    return MonthlyGridSeries(
        source_id=str(config["source_id"]),
        months=months,
        lat=lat,
        lon=lon,
        ewh_mm=field,
        valid_month=valid_month,
        month_status=month_status,
        metadata=metadata,
    )


def load_mascon(config: Mapping[str, Any]) -> MonthlyGridSeries:
    """Load an observed or reconstructed Mascon NetCDF into the common grid type."""
    return _load_grid(config, default_status="observed")


def load_custom_l3(config: Mapping[str, Any]) -> MonthlyGridSeries:
    """Load a local Level-3 NetCDF while honoring its monthly validity flag."""
    return _load_grid(config, default_status="observed")
