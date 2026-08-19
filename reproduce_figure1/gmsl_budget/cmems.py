from __future__ import annotations

from pathlib import Path

import numpy as np
import xarray as xr

from .masks import cell_area_weights
from .models import MonthlySeries, SpatialMask


def build_altimetry_mask(sla: xr.DataArray, ocean_fraction: xr.DataArray) -> SpatialMask:
    required = {"time", "latitude", "longitude"}
    if not required.issubset(sla.dims):
        raise ValueError("SLA must have time, latitude, and longitude dimensions")
    fraction = ocean_fraction.transpose("latitude", "longitude")
    if not np.array_equal(sla.latitude.values, fraction.latitude.values):
        raise ValueError("ocean fraction latitude does not match SLA")
    if not np.array_equal(sla.longitude.values, fraction.longitude.values):
        raise ValueError("ocean fraction longitude does not match SLA")
    values = sla.transpose("time", "latitude", "longitude").values
    support = np.all(np.isfinite(values), axis=0) & (fraction.values > 0.0)
    return SpatialMask(
        latitude=sla.latitude.values,
        longitude=sla.longitude.values,
        ocean_fraction=fraction.values,
        support=support,
        metadata={
            "name": "altimetry_global",
            "support_rule": "finite SLA in every input month and positive ocean fraction",
            "dynamic_monthly_mask": False,
        },
    )


def _find_sla(dataset: xr.Dataset, variable: str | None) -> xr.DataArray:
    if variable is not None:
        if variable not in dataset:
            raise ValueError(f"CMEMS variable not found: {variable}")
        return dataset[variable]
    if "sla" in dataset:
        return dataset["sla"]
    raise ValueError("CMEMS dataset does not contain an sla variable")


def _normalize_coordinates(sla: xr.DataArray) -> xr.DataArray:
    rename = {}
    if "lat" in sla.dims and "latitude" not in sla.dims:
        rename["lat"] = "latitude"
    if "lon" in sla.dims and "longitude" not in sla.dims:
        rename["lon"] = "longitude"
    sla = sla.rename(rename)
    if not {"time", "latitude", "longitude"}.issubset(sla.dims):
        raise ValueError("CMEMS SLA dimensions are not recognizable")
    if len(np.unique(sla.latitude.values)) != sla.sizes["latitude"]:
        raise ValueError("CMEMS latitude contains duplicates")
    normalized_longitude = np.mod(np.asarray(sla.longitude.values, dtype=np.float64), 360.0)
    if len(np.unique(np.round(normalized_longitude, 12))) != len(normalized_longitude):
        raise ValueError("CMEMS longitude contains duplicates")
    sla = sla.assign_coords(longitude=normalized_longitude)
    return sla.sortby("latitude").sortby("longitude")


def _unit_multiplier_to_mm(units: str) -> float:
    normalized = units.strip().lower()
    if normalized in {"m", "meter", "metre", "meters", "metres"}:
        return 1000.0
    if normalized in {"mm", "millimeter", "millimetre", "millimeters", "millimetres"}:
        return 1.0
    raise ValueError(f"unsupported SLA units: {units!r}")


def read_cmems_gmsl(
    path: str | Path,
    mask: SpatialMask,
    min_valid_weight_fraction: float,
    variable: str | None = None,
) -> MonthlySeries:
    if not 0 < min_valid_weight_fraction <= 1:
        raise ValueError("min_valid_weight_fraction must be in (0, 1]")
    source_path = Path(path)
    with xr.open_dataset(source_path) as dataset:
        sla = _normalize_coordinates(_find_sla(dataset, variable)).transpose(
            "time", "latitude", "longitude"
        )
        if not np.allclose(sla.latitude.values, mask.latitude):
            raise ValueError("CMEMS latitude does not match fixed mask")
        if not np.allclose(sla.longitude.values, np.mod(mask.longitude, 360.0)):
            raise ValueError("CMEMS longitude does not match fixed mask")
        multiplier = _unit_multiplier_to_mm(str(sla.attrs.get("units", "")))
        weights = cell_area_weights(mask.latitude, mask.longitude) * mask.ocean_fraction * mask.support
        total_weight = float(np.sum(weights))
        if not np.isfinite(total_weight) or total_weight <= 0:
            raise ValueError("fixed mask has no positive ocean weight")
        output = np.full(sla.sizes["time"], np.nan, dtype=np.float64)
        valid_fractions = np.zeros(sla.sizes["time"], dtype=np.float64)
        for index in range(sla.sizes["time"]):
            field = np.asarray(sla.isel(time=index).values, dtype=np.float64) * multiplier
            valid = np.isfinite(field) & mask.support
            valid_weight = float(np.sum(weights[valid]))
            valid_fractions[index] = valid_weight / total_weight
            if valid_fractions[index] >= min_valid_weight_fraction:
                output[index] = float(np.sum(np.where(valid, field, 0.0) * weights) / total_weight)
        time = sla.time.values
        attrs = dict(dataset.attrs)
        variable_name = str(sla.name)
    return MonthlySeries(
        time=time,
        values=output,
        name="gmsl_raw",
        units="mm",
        metadata={
            "source_path": str(source_path.resolve()),
            "source_variable": variable_name,
            "source_global_attributes": attrs,
            "input_units": str(sla.attrs.get("units", "")),
            "fixed_mask_name": mask.metadata.get("name"),
            "min_valid_weight_fraction": float(min_valid_weight_fraction),
            "valid_weight_fraction": valid_fractions.tolist(),
            "dynamic_monthly_mask": False,
            "gia_corrected": False,
        },
    )
