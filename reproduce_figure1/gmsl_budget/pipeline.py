from __future__ import annotations

from dataclasses import asdict, dataclass
from datetime import datetime, timezone
from hashlib import sha256
import json
from pathlib import Path
import shutil
from typing import Any, Mapping
from uuid import uuid4

import numpy as np
import pandas as pd
import xarray as xr

from .cmems import build_altimetry_mask, read_cmems_gmsl
from .config import PipelineConfig
from .gia import apply_scalar_gia
from .icgem import parse_gfc
from .masks import buffer_ocean_mask, load_cdt_coast_distance, load_cdt_ocean_fraction
from .models import MonthlySeries, SpatialMask, TrendResult
from .obd import area_average_mass_and_obd, preprocess_grace_coefficients
from .provenance import sha256_file, sha256_mask, software_versions, write_json
from .report import write_diagnostics, write_run_report
from .trend import fit_trend


@dataclass(frozen=True)
class PipelineRun:
    series: Mapping[str, MonthlySeries]
    masks: Mapping[str, SpatialMask]
    provenance: Mapping[str, Any]


def _input_inventory(config: PipelineConfig) -> tuple[dict[str, str], str]:
    low_degree = config.sagea_root / "data" / "L2_low_degrees"
    candidates = [
        config.cmems_sla_path,
        config.cdt_land_mask_path,
        config.cdt_distance_path,
        low_degree / "TN-11_C20_SLR_RL06.txt",
        low_degree / "TN-13_GEOC_CSR_RL06.2.txt",
        low_degree / "TN-14_C30_C20_SLR_GSFC.txt",
        config.sagea_root / "data" / "GIA" / "GIA.Caron_et_al_2018.txt",
        config.sagea_root / "data" / "auxiliary" / "LoveNumber.mat",
        *sorted(config.grace_gfc_dir.glob("*.gfc")),
    ]
    for optional in (config.cmems_indicator_path, config.steric_path):
        if optional is not None:
            candidates.append(optional)
    inventory = {
        str(path.resolve()): sha256_file(path) if path.is_file() else "missing"
        for path in candidates
    }
    serialized = json.dumps(inventory, ensure_ascii=False, sort_keys=True, separators=(",", ":")).encode("utf-8")
    return inventory, sha256(serialized).hexdigest()


def _serializable_attribute(value: Any) -> str | int | float:
    if isinstance(value, (str, int, float, np.integer, np.floating)) and not isinstance(value, bool):
        return value
    if isinstance(value, bool):
        return int(value)
    return json.dumps(value, ensure_ascii=False, sort_keys=True, default=str)


def _mask_with_hash(mask: SpatialMask, name: str) -> SpatialMask:
    metadata = dict(mask.metadata)
    metadata.update({"name": name, "sha256": sha256_mask(mask)})
    return SpatialMask(mask.latitude, mask.longitude, mask.ocean_fraction, mask.support, metadata)


def _open_sla_for_mask(config: PipelineConfig) -> xr.DataArray:
    dataset = xr.open_dataset(config.cmems_sla_path)
    try:
        variable = "sla" if "sla" in dataset else None
        if variable is None:
            raise ValueError("CMEMS dataset does not contain an sla variable")
        sla = dataset[variable]
        rename = {}
        if "lat" in sla.dims:
            rename["lat"] = "latitude"
        if "lon" in sla.dims:
            rename["lon"] = "longitude"
        sla = sla.rename(rename).assign_coords(longitude=np.mod(sla.longitude, 360.0)).sortby("latitude").sortby("longitude")
        periods = pd.DatetimeIndex(sla.time.values).to_period("M")
        keep = (periods >= pd.Period(config.start_month, "M")) & (periods <= pd.Period(config.end_month, "M"))
        if not np.any(keep):
            raise ValueError("CMEMS data do not overlap the configured date range")
        return sla.isel(time=np.flatnonzero(keep)).load()
    finally:
        dataset.close()


def _compute_run(config: PipelineConfig) -> PipelineRun:
    if config.steric_path is not None:
        raise NotImplementedError("qualified steric ingestion is intentionally deferred")
    sla = _open_sla_for_mask(config)
    ocean_fraction = load_cdt_ocean_fraction(
        config.cdt_land_mask_path,
        sla.latitude.values,
        sla.longitude.values,
        config.grid_spacing_degrees,
    )
    altimetry_mask = _mask_with_hash(build_altimetry_mask(sla, ocean_fraction), "altimetry_global")
    coast_distance = load_cdt_coast_distance(
        config.cdt_distance_path, sla.latitude.values, sla.longitude.values
    )
    budget_mask = _mask_with_hash(
        buffer_ocean_mask(altimetry_mask, coast_distance, config.coastal_buffer_km),
        "budget_common",
    )

    raw_global = read_cmems_gmsl(
        config.cmems_sla_path, altimetry_mask, config.min_valid_weight_fraction
    )
    raw_common_base = read_cmems_gmsl(
        config.cmems_sla_path, budget_mask, config.min_valid_weight_fraction
    )
    wanted = (raw_global.time.to_period("M") >= pd.Period(config.start_month, "M")) & (
        raw_global.time.to_period("M") <= pd.Period(config.end_month, "M")
    )
    raw_global = MonthlySeries(
        raw_global.time[wanted], raw_global.values[wanted], "gmsl_raw", "mm", {**raw_global.metadata, "mask_hash": altimetry_mask.metadata["sha256"]}
    )
    corrected_global = apply_scalar_gia(raw_global, config.altimetry_gia_rate_mm_per_year)
    common_wanted = (raw_common_base.time.to_period("M") >= pd.Period(config.start_month, "M")) & (
        raw_common_base.time.to_period("M") <= pd.Period(config.end_month, "M")
    )
    raw_common = MonthlySeries(
        raw_common_base.time[common_wanted],
        raw_common_base.values[common_wanted],
        "gmsl_budget_common_raw",
        "mm",
        {**raw_common_base.metadata, "mask_hash": budget_mask.metadata["sha256"]},
    )
    corrected_common_base = apply_scalar_gia(raw_common, config.altimetry_gia_rate_mm_per_year)
    corrected_common = MonthlySeries(
        corrected_common_base.time,
        corrected_common_base.values,
        "gmsl_budget_common_gia_corrected",
        "mm",
        corrected_common_base.metadata,
    )

    start = pd.Period(config.start_month, "M")
    end = pd.Period(config.end_month, "M")
    parsed = [parse_gfc(path, lmax=60) for path in sorted(config.grace_gfc_dir.glob("*.gfc"))]
    epochs = [item for item in parsed if start <= pd.Period(item.midpoint, "M") <= end]
    if not epochs:
        raise FileNotFoundError(f"no GFC epochs found in {config.grace_gfc_dir}")
    processed = preprocess_grace_coefficients(epochs, config.sagea_root, config.mass_preprocessing)
    mass, obd = area_average_mass_and_obd(
        processed,
        config.sagea_root,
        budget_mask,
        config.rho_freshwater_kg_m3,
        config.rho_seawater_kg_m3,
    )

    input_paths = [
        config.cmems_sla_path,
        config.cdt_land_mask_path,
        config.cdt_distance_path,
        *(item.path for item in epochs),
    ]
    provenance = {
        "input_sha256": {str(path.resolve()): sha256_file(path) for path in input_paths},
        "grace_preprocessing_hash": processed.preprocessing_hash,
        "mask_hashes": {
            "altimetry_global": altimetry_mask.metadata["sha256"],
            "budget_common": budget_mask.metadata["sha256"],
        },
    }
    return PipelineRun(
        series={
            item.name: item
            for item in (raw_global, corrected_global, raw_common, corrected_common, mass, obd)
        },
        masks={"altimetry_global": altimetry_mask, "budget_common": budget_mask},
        provenance=provenance,
    )


def _monthly_dataset(series: Mapping[str, MonthlySeries]) -> xr.Dataset:
    variables = {}
    for name, item in series.items():
        variables[name] = xr.DataArray(
            item.values,
            dims=("time",),
            coords={"time": item.time},
            attrs={"units": item.units, **{key: _serializable_attribute(value) for key, value in item.metadata.items()}},
        )
    return xr.Dataset(variables).sortby("time")


def _write_mask(mask: SpatialMask, path: Path) -> None:
    dataset = xr.Dataset(
        {
            "ocean_fraction": (("latitude", "longitude"), mask.ocean_fraction, {"units": "1"}),
            "support": (("latitude", "longitude"), mask.support.astype(np.int8), {"units": "1"}),
        },
        coords={"latitude": mask.latitude, "longitude": mask.longitude},
        attrs={key: _serializable_attribute(value) for key, value in mask.metadata.items()},
    )
    dataset.to_netcdf(path)


def _trend_rows(series: Mapping[str, MonthlySeries], hac_lags: int) -> tuple[list[TrendResult], list[str]]:
    results: list[TrendResult] = []
    warnings: list[str] = []
    for item in series.values():
        try:
            results.append(fit_trend(item, hac_lags))
        except ValueError as error:
            warnings.append(f"Trend unavailable for {item.name}: {error}")
    obd_trend = next((item.trend_mm_per_year for item in results if item.series_name == "obd"), None)
    if obd_trend is not None and not -0.30 <= obd_trend <= 0.10:
        warnings.append(f"OBD trend {obd_trend:.3f} mm/yr is outside the diagnostic range [-0.30, +0.10].")
    return results, warnings


def write_run_outputs(run: PipelineRun, output_dir: str | Path, config: PipelineConfig) -> None:
    destination = Path(output_dir)
    destination.mkdir(parents=True, exist_ok=False)
    dataset = _monthly_dataset(run.series)
    dataset.attrs.update(
        {
            "title": "Auditable full-ocean GMSL budget components",
            "configuration_sha256": config.sha256(),
            "budget_closure_available": 0,
        }
    )
    dataset.to_netcdf(destination / "monthly_budget.nc")
    dataset.to_dataframe().to_csv(destination / "monthly_budget.csv", index_label="time")
    trends, warnings = _trend_rows(run.series, config.hac_lags)
    pd.DataFrame([asdict(item) for item in trends]).to_csv(destination / "trend_summary.csv", index=False)
    _write_mask(run.masks["altimetry_global"], destination / "mask_altimetry_global.nc")
    _write_mask(run.masks["budget_common"], destination / "mask_budget_common.nc")
    write_json(destination / "config_resolved.json", config.resolved_dict())
    provenance = {
        **dict(run.provenance),
        "configuration_sha256": config.sha256(),
        "budget_closure_available": False,
        "steric_status": "not supplied; no closure computed",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "software_versions": software_versions(),
        "scientific_warnings": warnings,
    }
    write_json(destination / "provenance.json", provenance)
    write_diagnostics(run.series, destination / "diagnostics.png")
    write_run_report(
        destination / "run_report.md",
        destination.name,
        config.sha256(),
        trends,
        False,
        warnings,
    )


def run_pipeline(config: PipelineConfig) -> Path:
    config_hash = config.sha256()
    inventory, inventory_hash = _input_inventory(config)
    run_id = config.run_id or f"gmsl-{config_hash[:8]}-{inventory_hash[:8]}"
    output_root = config.output_root
    output_root.mkdir(parents=True, exist_ok=True)
    destination = output_root / run_id
    if destination.exists():
        provenance_path = destination / "provenance.json"
        if provenance_path.is_file():
            previous = json.loads(provenance_path.read_text(encoding="utf-8"))
            if (
                previous.get("configuration_sha256") == config_hash
                and previous.get("input_inventory_sha256") == inventory_hash
            ):
                return destination
            if previous.get("configuration_sha256") == config_hash:
                raise FileExistsError(f"run directory exists with a different input inventory hash: {destination}")
        raise FileExistsError(f"run directory exists with a different configuration hash: {destination}")
    computed = _compute_run(config)
    run = PipelineRun(
        computed.series,
        computed.masks,
        {
            **dict(computed.provenance),
            "input_inventory": inventory,
            "input_inventory_sha256": inventory_hash,
        },
    )
    staging = output_root / f".{run_id}.tmp-{uuid4().hex}"
    try:
        write_run_outputs(run, staging, config)
        staging.replace(destination)
    except Exception:
        if staging.exists():
            shutil.rmtree(staging)
        raise
    return destination
