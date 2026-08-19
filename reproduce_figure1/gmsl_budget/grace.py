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
        with xr.open_dataset(path) as source:
            dataset = source.load()
    except ValueError as error:
        if "backends" not in str(error):
            raise
        dataset = _read_hdf5_netcdf(path, variable)
    candidates = [variable] if variable else ["lwe_thickness", "ewh", "EWH"]
    variable_name = next((name for name in candidates if name and name in dataset), None)
    if variable_name is None:
        raise ValueError("mascon file does not contain a recognized EWH variable")
    return dataset, variable_name


def _verify_center(center: str, attrs: Mapping[str, object]) -> None:
    identity = " ".join(str(value) for value in attrs.values()).lower()
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
    field = dataset[variable_name].rename(rename).transpose("time", "latitude", "longitude")
    field = field.assign_coords(longitude=np.mod(field.longitude.values, 360.0)).sortby("latitude").sortby("longitude")
    if not np.allclose(field.latitude.values, mask.latitude):
        raise ValueError("mascon latitude does not match fixed mask")
    if not np.allclose(field.longitude.values, np.mod(mask.longitude, 360.0)):
        raise ValueError("mascon longitude does not match fixed mask")
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
    preprocessing_payload = {
        "center": center.upper(),
        "source": str(dataset.attrs.get("source", "")),
        "product_version": str(dataset.attrs.get("product_version", "")),
        "gia": str(dataset.attrs.get("GIA_removed", dataset.attrs.get("GIA_Removed", ""))),
    }
    preprocessing_hash = sha256(
        json.dumps(preprocessing_payload, sort_keys=True, separators=(",", ":")).encode("utf-8")
    ).hexdigest()
    return MonthlySeries(
        field.time.values,
        values,
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
            "dynamic_monthly_mask": False,
        },
    )
