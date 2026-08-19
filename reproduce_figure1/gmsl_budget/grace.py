from __future__ import annotations

from hashlib import sha256
import json
from pathlib import Path
import re
from typing import Mapping, Sequence

import h5py
import numpy as np
import pandas as pd
import xarray as xr

from .masks import cell_area_weights
from .models import MonthlySeries, SpatialMask


def ewh_to_sea_level(
    ewh: np.ndarray,
    rho_freshwater: float = 1000.0,
    rho_seawater: float = 1028.0,
) -> np.ndarray:
    if rho_freshwater <= 0 or rho_seawater <= 0:
        raise ValueError("water densities must be positive")
    return np.asarray(ewh, dtype=np.float64) * float(rho_freshwater) / float(rho_seawater)


def align_common_months(series: Sequence[MonthlySeries]) -> list[MonthlySeries]:
    if not series:
        raise ValueError("at least one monthly series is required")
    common = set(series[0].time.to_period("M"))
    for item in series[1:]:
        common &= set(item.time.to_period("M"))
    if not common:
        raise ValueError("monthly series have no common months")
    periods = pd.PeriodIndex(sorted(common), freq="M")
    output = []
    for item in series:
        lookup = {period: value for period, value in zip(item.time.to_period("M"), item.values)}
        values = np.asarray([lookup[period] for period in periods], dtype=np.float64)
        metadata = dict(item.metadata)
        metadata["aligned_to_common_months"] = True
        metadata["common_month_count"] = len(periods)
        output.append(MonthlySeries(periods.to_timestamp() + pd.Timedelta(days=14), values, item.name, item.units, metadata))
    return output


def ensemble_mean(series_by_center: Mapping[str, MonthlySeries]) -> MonthlySeries:
    if not series_by_center:
        raise ValueError("at least one GRACE center is required")
    centers = sorted(series_by_center)
    for center in centers:
        declared = str(series_by_center[center].metadata.get("center", "")).upper()
        if declared and declared != center.upper():
            raise ValueError(f"series center {declared} does not match mapping key {center}")
        if series_by_center[center].units != "mm":
            raise ValueError("all ensemble members must use mm")
    mask_hashes = {series_by_center[center].metadata.get("mask_hash") for center in centers}
    if len(mask_hashes) != 1:
        raise ValueError("ensemble members must share one mask_hash")
    aligned = align_common_months([series_by_center[center] for center in centers])
    matrix = np.vstack([item.values for item in aligned])
    values = np.mean(matrix, axis=0)
    if len(centers) > 1:
        spread = np.std(matrix, axis=0, ddof=1)
    else:
        spread = np.full(matrix.shape[1], np.nan)
    metadata = {
        "centers": centers,
        "mask_hash": next(iter(mask_hashes)),
        "preprocessing_hashes": {
            center: series_by_center[center].metadata.get("preprocessing_hash") for center in centers
        },
        "ensemble_spread_mm": spread.tolist(),
        "ensemble_method": "unweighted arithmetic mean on common months",
    }
    if all(series_by_center[center].metadata.get("gia_corrected") is True for center in centers):
        metadata.update(
            {
                "gia_corrected": True,
                "mascon_gia_policy": "use product correction; do not apply GIA again",
            }
        )
    return MonthlySeries(aligned[0].time, values, "ocean_mass_ensemble", "mm", metadata)


def _decode(value):
    if isinstance(value, bytes):
        return value.decode("utf-8", errors="replace")
    if isinstance(value, np.ndarray) and value.size == 1:
        return _decode(value.reshape(-1)[0])
    if isinstance(value, np.generic):
        return value.item()
    return value


def _read_hdf5_netcdf(path: Path, variable: str | None) -> xr.Dataset:
    with h5py.File(path, "r") as source:
        candidates = [variable] if variable else ["lwe_thickness", "ewh", "EWH"]
        variable_name = next((name for name in candidates if name and name in source), None)
        if variable_name is None:
            raise ValueError("mascon file does not contain a recognized EWH variable")
        lat_name = "lat" if "lat" in source else "latitude"
        lon_name = "lon" if "lon" in source else "longitude"
        if lat_name not in source or lon_name not in source or "time" not in source:
            raise ValueError("mascon file is missing time/latitude/longitude coordinates")
        latitude = np.asarray(source[lat_name], dtype=np.float64)
        longitude = np.asarray(source[lon_name], dtype=np.float64)
        time_values = np.asarray(source["time"], dtype=np.float64)
        time_attrs = {str(key): _decode(value) for key, value in source["time"].attrs.items()}
        time_units = str(time_attrs.get("units", time_attrs.get("Units", "")))
        match = re.match(r"days since ([0-9]{4}-[0-9]{2}-[0-9]{2})", time_units)
        if not match:
            raise ValueError(f"unsupported mascon time units: {time_units!r}")
        time = pd.Timestamp(match.group(1)) + pd.to_timedelta(time_values, unit="D")
        data = np.asarray(source[variable_name], dtype=np.float64)
        variable_attrs = {str(key): _decode(value) for key, value in source[variable_name].attrs.items()}
        fill = variable_attrs.get("_FillValue")
        if fill is not None:
            data[np.isclose(data, float(fill))] = np.nan
        global_attrs = {str(key): _decode(value) for key, value in source.attrs.items()}
    return xr.Dataset(
        {
            variable_name: (
                ("time", "latitude", "longitude"),
                data,
                {key: value for key, value in variable_attrs.items() if np.isscalar(value)},
            )
        },
        coords={"time": time, "latitude": latitude, "longitude": longitude},
        attrs={key: value for key, value in global_attrs.items() if np.isscalar(value)},
    )


def _open_mascon(path: Path, variable: str | None) -> tuple[xr.Dataset, str]:
    try:
        dataset = xr.open_dataset(path, decode_times=False)
    except ValueError as error:
        if "backends" not in str(error):
            raise
        dataset = _read_hdf5_netcdf(path, variable)
    candidates = [variable] if variable else ["lwe_thickness", "ewh", "EWH"]
    variable_name = next((name for name in candidates if name and name in dataset), None)
    if variable_name is None:
        raise ValueError("mascon file does not contain a recognized EWH variable")
    return dataset, variable_name


def _decode_mascon_time(dataset: xr.Dataset) -> pd.DatetimeIndex:
    coordinate = dataset["time"]
    if np.issubdtype(coordinate.dtype, np.datetime64):
        return pd.DatetimeIndex(pd.to_datetime(coordinate.values))
    units = str(coordinate.attrs.get("units", coordinate.attrs.get("Units", "")))
    match = re.match(r"days since ([0-9]{4}-[0-9]{2}-[0-9]{2})", units)
    if not match:
        raise ValueError(f"unsupported mascon time units: {units!r}")
    return pd.DatetimeIndex(
        pd.Timestamp(match.group(1)) + pd.to_timedelta(np.asarray(coordinate.values, dtype=np.float64), unit="D")
    )


def _collapse_monthly(time: pd.DatetimeIndex, values: np.ndarray) -> tuple[pd.DatetimeIndex, np.ndarray]:
    frame = pd.DataFrame({"month": time.to_period("M"), "value": np.asarray(values, dtype=np.float64)})
    monthly = frame.groupby("month", sort=True, observed=True)["value"].mean()
    return monthly.index.to_timestamp() + pd.Timedelta(days=14), monthly.to_numpy(dtype=np.float64)


def _reject_all_zero(values: np.ndarray, center: str) -> None:
    finite = np.asarray(values, dtype=np.float64)[np.isfinite(values)]
    if finite.size and float(np.max(np.abs(finite))) <= 1.0e-12:
        raise ValueError(f"{center} mascon scientific period is all-zero")


def _verify_center(center: str, attrs: Mapping[str, object]) -> None:
    identity_fields = (
        "title",
        "subtitle",
        "summary",
        "institution",
        "creator_institution",
        "publisher_institution",
        "source",
        "filename",
    )
    identity = " ".join(str(attrs.get(key, "")) for key in identity_fields).lower()
    aliases = {
        "CSR": ("csr", "center for space research"),
        "JPL": ("jpl", "jet propulsion laboratory"),
        "GFZ": ("gfz", "geoforschungszentrum"),
        "GSFC": ("gsfc", "goddard"),
    }
    expected = center.upper()
    if expected not in aliases:
        raise ValueError(f"unsupported GRACE center: {center}")
    detected = {name for name, tokens in aliases.items() if any(token in identity for token in tokens)}
    if expected not in detected or any(name != expected for name in detected):
        raise ValueError(f"declared center {expected} does not match file metadata: {sorted(detected)}")


def _ewh_units_to_mm(units: str) -> float:
    normalized = units.strip().lower()
    if normalized in {"mm", "millimeter", "millimetre"}:
        return 1.0
    if normalized in {"cm", "centimeter", "centimetre"}:
        return 10.0
    if normalized in {"m", "meter", "metre"}:
        return 1000.0
    raise ValueError(f"unsupported mascon EWH units: {units!r}")


def read_mascon_ocean_series(
    path: str | Path,
    center: str,
    mask: SpatialMask,
    variable: str | None = None,
    rho_freshwater: float = 1000.0,
    rho_seawater: float = 1028.0,
) -> MonthlySeries:
    source_path = Path(path)
    dataset, variable_name = _open_mascon(source_path, variable)
    _verify_center(center, dataset.attrs)
    rename = {}
    for old, new in (("lat", "latitude"), ("lon", "longitude")):
        if old in dataset[variable_name].dims and new not in dataset[variable_name].dims:
            rename[old] = new
    decoded_time = _decode_mascon_time(dataset)
    field = dataset[variable_name].rename(rename).transpose("time", "latitude", "longitude")
    field = field.assign_coords(
        time=decoded_time,
        longitude=np.mod(field.longitude.values, 360.0),
    ).sortby("latitude").sortby("longitude")
    field = field.sel(
        latitude=xr.DataArray(mask.latitude, dims="latitude"),
        longitude=xr.DataArray(np.mod(mask.longitude, 360.0), dims="longitude"),
        method="nearest",
    ).load()
    units = str(field.attrs.get("units", field.attrs.get("Units", "")))
    multiplier = _ewh_units_to_mm(units)
    weights = cell_area_weights(mask.latitude, mask.longitude) * mask.ocean_fraction * mask.support
    total_weight = float(np.sum(weights))
    if total_weight <= 0:
        raise ValueError("fixed mask has no positive ocean weight")
    values = np.full(field.sizes["time"], np.nan, dtype=np.float64)
    for index in range(field.sizes["time"]):
        ewh_mm = np.asarray(field.isel(time=index).values, dtype=np.float64) * multiplier
        valid = np.isfinite(ewh_mm) & mask.support
        if float(np.sum(weights[valid])) / total_weight >= 0.995:
            mean_ewh = float(np.sum(np.where(valid, ewh_mm, 0.0) * weights) / total_weight)
            values[index] = ewh_to_sea_level(mean_ewh, rho_freshwater, rho_seawater)
    monthly_time, monthly_values = _collapse_monthly(pd.DatetimeIndex(field.time.values), values)
    _reject_all_zero(monthly_values, center.upper())
    preprocessing_payload = {
        "center": center.upper(),
        "source": str(dataset.attrs.get("source", "")),
        "product_version": str(dataset.attrs.get("product_version", "")),
        "gia": str(dataset.attrs.get("GIA_removed", dataset.attrs.get("GIA_Removed", ""))),
    }
    preprocessing_hash = sha256(
        json.dumps(preprocessing_payload, sort_keys=True, separators=(",", ":")).encode("utf-8")
    ).hexdigest()
    dataset.close()
    return MonthlySeries(
        monthly_time,
        monthly_values,
        f"ocean_mass_{center.lower()}",
        "mm",
        {
            "center": center.upper(),
            "source_path": str(source_path.resolve()),
            "source_variable": variable_name,
            "source_units": units,
            "rho_freshwater_kg_m3": float(rho_freshwater),
            "rho_seawater_kg_m3": float(rho_seawater),
            "mask_hash": mask.metadata.get("sha256"),
            "preprocessing_hash": preprocessing_hash,
            "preprocessing": preprocessing_payload,
            "gia_corrected": True,
            "mascon_gia_policy": "use product correction; do not apply GIA again",
            "dynamic_monthly_mask": False,
            "spatial_interpolation": "nearest neighbour to fixed mask grid",
            "duplicate_month_policy": "arithmetic mean",
        },
    )


def _hdf5_vector(group: h5py.File, key: str) -> np.ndarray:
    return np.asarray(group[key], dtype=np.float64).reshape(-1)


def read_gsfc_sla_mascon_series(
    path: str | Path,
    mask: SpatialMask,
    rho_freshwater: float = 1000.0,
    rho_seawater: float = 1028.0,
) -> MonthlySeries:
    source_path = Path(path)
    with h5py.File(source_path, "r") as source:
        calendar = np.asarray(source["/time/yyyy_doy_yrplot_middle"], dtype=np.float64)
        if calendar.shape[0] != 3 and calendar.shape[1] == 3:
            calendar = calendar.T
        if calendar.shape[0] != 3:
            raise ValueError("GSFC HDF5 time table must have year, day-of-year, and decimal year")
        years = calendar[0].astype(int)
        days = calendar[1].astype(int)
        time = pd.DatetimeIndex(
            [pd.Timestamp(year=int(year), month=1, day=1) + pd.Timedelta(days=int(day) - 1) for year, day in zip(years, days)]
        )
        area = _hdf5_vector(source, "/mascon/area_km2")
        location = _hdf5_vector(source, "/mascon/location")
        basin = _hdf5_vector(source, "/mascon/basin")
        latitude = _hdf5_vector(source, "/mascon/lat_center")
        longitude = np.mod(_hdf5_vector(source, "/mascon/lon_center"), 360.0)
        cmwe = np.asarray(source["/solution/cmwe"], dtype=np.float64)
    if cmwe.shape == (len(time), len(area)):
        cmwe = cmwe.T
    if cmwe.shape != (len(area), len(time)):
        raise ValueError("GSFC HDF5 solution dimensions do not match mascon and time tables")
    lat_index = np.asarray([int(np.argmin(np.abs(mask.latitude - value))) for value in latitude], dtype=int)
    lon_index = np.asarray(
        [int(np.argmin(np.abs((np.mod(mask.longitude, 360.0) - value + 180.0) % 360.0 - 180.0))) for value in longitude],
        dtype=int,
    )
    in_fixed_mask = mask.support[lat_index, lon_index] & (mask.ocean_fraction[lat_index, lon_index] > 0.0)
    selected = (location == 90.0) & (basin == 0.0) & in_fixed_mask & np.isfinite(area) & (area > 0.0)
    if not np.any(selected):
        raise ValueError("GSFC HDF5 has no open-ocean mascons inside the fixed mask")
    mean_cmwe = np.sum(cmwe[selected] * area[selected, None], axis=0) / float(np.sum(area[selected]))
    values = ewh_to_sea_level(mean_cmwe * 10.0, rho_freshwater, rho_seawater)
    monthly_time, monthly_values = _collapse_monthly(time, values)
    _reject_all_zero(monthly_values, "GSFC")
    preprocessing_payload = {
        "center": "GSFC",
        "product": "RL06v2.0 SLA-ICE6GD native equal-area mascons",
        "gia": "ICE6G-D removed",
        "gad": "restored with global ocean mean removed",
    }
    preprocessing_hash = sha256(
        json.dumps(preprocessing_payload, sort_keys=True, separators=(",", ":")).encode("utf-8")
    ).hexdigest()
    return MonthlySeries(
        monthly_time,
        monthly_values,
        "ocean_mass_gsfc",
        "mm",
        {
            "center": "GSFC",
            "source_path": str(source_path.resolve()),
            "source_variable": "/solution/cmwe",
            "source_units": "cm equivalent water height",
            "rho_freshwater_kg_m3": float(rho_freshwater),
            "rho_seawater_kg_m3": float(rho_seawater),
            "mask_hash": mask.metadata.get("sha256"),
            "preprocessing_hash": preprocessing_hash,
            "preprocessing": preprocessing_payload,
            "gia_corrected": True,
            "mascon_gia_policy": "use product correction; do not apply GIA again",
            "native_ocean_mascon_count": int(np.sum(selected)),
            "spatial_interpolation": "fixed mask sampled at native mascon centres",
            "duplicate_month_policy": "arithmetic mean",
        },
    )
