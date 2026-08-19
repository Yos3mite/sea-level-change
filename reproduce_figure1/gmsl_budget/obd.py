from __future__ import annotations

from contextlib import contextmanager
from dataclasses import dataclass
from datetime import date
from hashlib import sha256
import json
import os
from pathlib import Path
import sys
from typing import Any, Mapping, Sequence

import numpy as np
import pandas as pd
from scipy.io import loadmat
import xarray as xr

from .icgem import GfcEpoch
from .masks import cell_area_weights
from .models import MonthlySeries, SpatialMask


EARTH_RADIUS_M = 6_378_136.3
EARTH_DENSITY_KG_M3 = 5_517.0
FRESHWATER_DENSITY_KG_M3 = 1_000.0


def _readonly(values: np.ndarray) -> np.ndarray:
    array = np.array(values, dtype=np.float64, copy=True)
    array.setflags(write=False)
    return array


@dataclass(frozen=True)
class ProcessedCoefficients:
    c: np.ndarray
    s: np.ndarray
    start: tuple[date, ...]
    end: tuple[date, ...]
    midpoint: tuple[date, ...]
    lmax: int
    preprocessing_hash: str
    metadata: Mapping[str, Any]

    def __post_init__(self) -> None:
        c = _readonly(self.c)
        s = _readonly(self.s)
        expected = (len(self.midpoint), self.lmax + 1, self.lmax + 1)
        if c.shape != expected or s.shape != expected:
            raise ValueError(f"coefficient arrays must have shape {expected}")
        if len(self.start) != len(self.midpoint) or len(self.end) != len(self.midpoint):
            raise ValueError("coefficient dates must have equal lengths")
        months = pd.PeriodIndex(self.midpoint, freq="M")
        if months.duplicated().any():
            raise ValueError("duplicate month in processed coefficients")
        if not self.preprocessing_hash:
            raise ValueError("preprocessing_hash must not be empty")
        object.__setattr__(self, "c", c)
        object.__setattr__(self, "s", s)
        object.__setattr__(self, "metadata", dict(self.metadata))


def load_wang_love_numbers(sagea_root: str | Path, lmax: int) -> tuple[np.ndarray, np.ndarray]:
    if lmax < 0 or lmax > 360:
        raise ValueError("Wang Love numbers support 0 <= lmax <= 360")
    path = Path(sagea_root) / "data" / "auxiliary" / "LoveNumber.mat"
    love = np.asarray(loadmat(path)["love"], dtype=np.float64)
    if love.shape[0] < lmax or love.shape[1] < 4:
        raise ValueError("LoveNumber.mat does not cover the requested degree")
    h = np.concatenate(([0.0], love[:lmax, 1]))
    k = np.concatenate(([0.0], love[:lmax, 3]))
    return k.astype(np.float64), h.astype(np.float64)


def vertical_conversion_weights(
    lmax: int,
    k: np.ndarray,
    h: np.ndarray,
    radius_m: float = EARTH_RADIUS_M,
) -> np.ndarray:
    k = np.asarray(k, dtype=np.float64)
    h = np.asarray(h, dtype=np.float64)
    if k.shape[0] < lmax + 1 or h.shape[0] < lmax + 1:
        raise ValueError("Love-number arrays do not cover lmax")
    denominator = 1.0 + k[: lmax + 1]
    if np.any(np.isclose(denominator, 0.0)):
        raise ValueError("vertical conversion is singular because 1 + k is zero")
    return float(radius_m) * h[: lmax + 1] / denominator


def ewh_conversion_weights(
    lmax: int,
    k: np.ndarray,
    radius_m: float = EARTH_RADIUS_M,
    earth_density: float = EARTH_DENSITY_KG_M3,
    water_density: float = FRESHWATER_DENSITY_KG_M3,
) -> np.ndarray:
    k = np.asarray(k, dtype=np.float64)
    if k.shape[0] < lmax + 1:
        raise ValueError("Love-number array does not cover lmax")
    if earth_density <= 0 or water_density <= 0:
        raise ValueError("densities must be positive")
    degree = np.arange(lmax + 1, dtype=np.float64)
    denominator = 1.0 + k[: lmax + 1]
    if np.any(np.isclose(denominator, 0.0)):
        raise ValueError("EWH conversion is singular because 1 + k is zero")
    return (2.0 * degree + 1.0) / denominator * float(radius_m) * float(earth_density) / (3.0 * float(water_density))


@contextmanager
def _sagea_context(sagea_root: Path):
    root = str(sagea_root.resolve())
    old_cwd = os.getcwd()
    inserted = root not in sys.path
    if inserted:
        sys.path.insert(0, root)
    os.chdir(root)
    try:
        yield
    finally:
        os.chdir(old_cwd)
        if inserted:
            sys.path.remove(root)


def _file_sha256(path: Path) -> str:
    digest = sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def preprocess_grace_coefficients(
    epochs: Sequence[GfcEpoch],
    sagea_root: str | Path,
    preprocessing: Mapping[str, Any],
) -> ProcessedCoefficients:
    if not epochs:
        raise ValueError("at least one GFC epoch is required")
    ordered = sorted(epochs, key=lambda item: item.midpoint)
    lmax_values = {item.lmax for item in ordered}
    if len(lmax_values) != 1:
        raise ValueError("all GFC epochs must share lmax")
    lmax = next(iter(lmax_values))
    root = Path(sagea_root)
    low_degree_dir = root / "data" / "L2_low_degrees"
    low_degree_paths = (
        low_degree_dir / "TN-11_C20_SLR_RL06.txt",
        low_degree_dir / "TN-13_GEOC_CSR_RL06.2.txt",
        low_degree_dir / "TN-14_C30_C20_SLR_GSFC.txt",
    )
    gia_path = root / "data" / "GIA" / "GIA.Caron_et_al_2018.txt"
    required = (*low_degree_paths, gia_path, root / "data" / "auxiliary" / "LoveNumber.mat")
    missing = [str(path) for path in required if not path.is_file()]
    if missing:
        raise FileNotFoundError(f"SaGEA preprocessing inputs are missing: {missing}")
    c = np.stack([item.c for item in ordered]).astype(np.float64)
    s = np.stack([item.s for item in ordered]).astype(np.float64)
    starts = tuple(item.start for item in ordered)
    ends = tuple(item.end for item in ordered)
    midpoints = tuple(item.midpoint for item in ordered)
    with _sagea_context(root):
        from pysrc.auxiliary.load_file.LoadL2LowDeg import load_low_degs
        from pysrc.auxiliary.load_file.LoadL2SH import load_SHC
        from pysrc.data_class.SHC import SHC

        shc = SHC(c, s)
        shc.dates = [list(starts), list(ends)]
        low_degrees = load_low_degs(*low_degree_paths)
        shc.replace_low_degs(list(starts), list(ends), low_deg=low_degrees, deg1=True, c20=True, c30=True)
        gia_trend = load_SHC(gia_path, key="", lmax=lmax)
        shc.subtract(gia_trend.expand(list(midpoints)))
        shc.de_background()
        processed_c, processed_s = shc.get_cs2d()
    input_hashes = {str(item.path): _file_sha256(item.path) for item in ordered}
    auxiliary_hashes = {str(path): _file_sha256(path) for path in required}
    payload = {
        "preprocessing": dict(preprocessing),
        "input_hashes": input_hashes,
        "auxiliary_hashes": auxiliary_hashes,
        "lmax": lmax,
        "months": [str(pd.Period(value, freq="M")) for value in midpoints],
        "operations": ["replace degree-1/C20/C30", "subtract Caron2018 GIA", "remove temporal mean"],
        "filter": "ICGEM DDK1 retained; no second Gaussian filter",
    }
    preprocessing_hash = sha256(
        json.dumps(payload, ensure_ascii=False, sort_keys=True, separators=(",", ":")).encode("utf-8")
    ).hexdigest()
    return ProcessedCoefficients(
        processed_c,
        processed_s,
        starts,
        ends,
        midpoints,
        lmax,
        preprocessing_hash,
        payload,
    )


def _synthesize(
    processed: ProcessedCoefficients,
    sagea_root: str | Path,
    latitude: np.ndarray,
    longitude: np.ndarray,
    degree_weights: np.ndarray,
    name: str,
    units: str,
) -> xr.DataArray:
    latitude = np.asarray(latitude, dtype=np.float64)
    longitude = np.asarray(longitude, dtype=np.float64)
    output = np.empty((len(processed.midpoint), len(latitude), len(longitude)), dtype=np.float64)
    with _sagea_context(Path(sagea_root)):
        from pysrc.post_processing.harmonic.Harmonic import Harmonic

        harmonic = Harmonic(latitude, longitude, processed.lmax, option=1)
        for index in range(len(processed.midpoint)):
            weighted_c = processed.c[index] * degree_weights[:, None]
            weighted_s = processed.s[index] * degree_weights[:, None]
            output[index] = harmonic.synthesis(weighted_c, weighted_s)
    return xr.DataArray(
        output,
        dims=("time", "latitude", "longitude"),
        coords={
            "time": pd.DatetimeIndex(processed.midpoint),
            "latitude": latitude,
            "longitude": longitude,
        },
        name=name,
        attrs={
            "units": units,
            "preprocessing_hash": processed.preprocessing_hash,
            "harmonic_synthesis": "SaGEA Harmonic, fully normalized coefficients",
        },
    )


def convert_coefficients_to_vertical_displacement(
    processed: ProcessedCoefficients,
    sagea_root: str | Path,
    latitude: np.ndarray,
    longitude: np.ndarray,
) -> xr.DataArray:
    k, h = load_wang_love_numbers(sagea_root, processed.lmax)
    weights = vertical_conversion_weights(processed.lmax, k, h)
    result = _synthesize(processed, sagea_root, latitude, longitude, weights, "vertical_displacement", "m")
    result.attrs["sign_convention"] = "upward positive; subsidence negative"
    result.attrs["conversion"] = "R * h_l / (1 + k_l)"
    return result


def convert_coefficients_to_ewh(
    processed: ProcessedCoefficients,
    sagea_root: str | Path,
    latitude: np.ndarray,
    longitude: np.ndarray,
) -> xr.DataArray:
    k, _ = load_wang_love_numbers(sagea_root, processed.lmax)
    weights = ewh_conversion_weights(processed.lmax, k)
    return _synthesize(processed, sagea_root, latitude, longitude, weights, "equivalent_water_height", "m")


def assert_shared_preprocessing(mass_hash: str, obd_hash: str) -> None:
    if mass_hash != obd_hash:
        raise ValueError("mass and OBD preprocessing hash must match")


def area_average_obd(
    displacement: xr.DataArray,
    mask: SpatialMask,
    preprocessing_hash: str,
) -> MonthlySeries:
    if str(displacement.attrs.get("units", "")) != "m":
        raise ValueError("vertical displacement must use metres")
    assert_shared_preprocessing(preprocessing_hash, str(displacement.attrs.get("preprocessing_hash", "")))
    grid = displacement.transpose("time", "latitude", "longitude")
    if not np.allclose(grid.latitude.values, mask.latitude) or not np.allclose(grid.longitude.values, mask.longitude):
        raise ValueError("vertical displacement grid does not match fixed mask")
    weights = cell_area_weights(mask.latitude, mask.longitude) * mask.ocean_fraction * mask.support
    total_weight = float(np.sum(weights))
    if total_weight <= 0:
        raise ValueError("fixed mask has no positive ocean weight")
    values = np.full(grid.sizes["time"], np.nan, dtype=np.float64)
    for index in range(grid.sizes["time"]):
        field = np.asarray(grid.isel(time=index).values, dtype=np.float64)
        valid = np.isfinite(field) & mask.support
        if float(np.sum(weights[valid])) / total_weight >= 0.995:
            values[index] = float(np.sum(np.where(valid, field, 0.0) * weights) / total_weight * 1000.0)
    return MonthlySeries(
        grid.time.values,
        values,
        "obd",
        "mm",
        {
            "preprocessing_hash": preprocessing_hash,
            "mask_hash": mask.metadata.get("sha256"),
            "sign_convention": "upward positive; subsidence negative",
            "source_quantity": "GRACE/GRACE-FO elastic vertical displacement",
        },
    )
