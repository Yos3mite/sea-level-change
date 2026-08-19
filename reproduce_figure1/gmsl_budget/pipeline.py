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
from .figure1_mass import build_ssa_ensemble
from .gia import apply_piecewise_trend_correction, apply_scalar_gia
from .grace import read_gsfc_sla_mascon_series, read_mascon_ocean_series
from .icgem import parse_gfc
from .masks import (
    buffer_ocean_mask,
    load_cdt_coast_distance,
    load_cdt_ocean_fraction,
    resample_mask_nearest,
)
from .models import MonthlySeries, SpatialMask, TrendResult
from .obd import area_average_mass_and_obd, preprocess_grace_coefficients
from .provenance import sha256_file, sha256_mask, software_versions, write_json
from .report import write_diagnostics, write_mass_figure, write_run_report
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
    for optional in (
        config.cmems_indicator_path,
        config.steric_path,
        config.csr_mascon_path,
        config.jpl_mascon_path,
        config.gsfc_sla_mascon_path,
    ):
        if optional is not None:
            candidates.append(optional)
    inventory = {
        str(path.resolve()): sha256_file(path) if path.is_file() else "missing"
        for path in candidates
    }
    serialized = json.dumps(inventory, ensure_ascii=False, sort_keys=True, separators=(",", ":")).encode("utf-8")
    return inventory, sha256(serialized).hexdigest()


def _source_inventory() -> tuple[dict[str, str], str]:
    package_root = Path(__file__).resolve().parent
    paths = [*sorted(package_root.glob("*.py")), package_root.parent / "run_gmsl_budget.py"]
    inventory = {str(path.resolve()): sha256_file(path) for path in paths}
    serialized = json.dumps(inventory, ensure_ascii=False, sort_keys=True, separators=(",", ":")).encode("utf-8")
    return inventory, sha256(serialized).hexdigest()


def _default_run_id(config_hash: str, inventory_hash: str, source_hash: str) -> str:
    return f"gmsl-{config_hash[:8]}-{inventory_hash[:8]}-{source_hash[:8]}"


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


def _subset_months(series: MonthlySeries, start_month: str, end_month: str, name: str | None = None) -> MonthlySeries:
    periods = series.time.to_period("M")
    keep = (periods >= pd.Period(start_month, "M")) & (periods <= pd.Period(end_month, "M"))
    if not np.any(keep):
        raise ValueError(f"{series.name} does not overlap {start_month} through {end_month}")
    return MonthlySeries(
        series.time[keep],
        series.values[keep],
        name or series.name,
        series.units,
        series.metadata,
    )


def _renamed(series: MonthlySeries, name: str, **metadata: Any) -> MonthlySeries:
    return MonthlySeries(series.time, series.values, name, series.units, {**series.metadata, **metadata})


def _read_mascon_centers(config: PipelineConfig, mask: SpatialMask) -> dict[str, MonthlySeries]:
    paths = (config.csr_mascon_path, config.jpl_mascon_path, config.gsfc_sla_mascon_path)
    if any(path is None for path in paths):
        raise ValueError("CSR, JPL, and GSFC Mascon paths must all be supplied together")
    return {
        "CSR": _subset_months(
            read_mascon_ocean_series(
                config.csr_mascon_path, "CSR", mask,
                rho_freshwater=config.rho_freshwater_kg_m3,
                rho_seawater=config.rho_seawater_kg_m3,
            ),
            config.start_month,
            config.end_month,
        ),
        "JPL": _subset_months(
            read_mascon_ocean_series(
                config.jpl_mascon_path, "JPL", mask,
                rho_freshwater=config.rho_freshwater_kg_m3,
                rho_seawater=config.rho_seawater_kg_m3,
            ),
            config.start_month,
            config.end_month,
        ),
        "GSFC": _subset_months(
            read_gsfc_sla_mascon_series(
                config.gsfc_sla_mascon_path, mask,
                rho_freshwater=config.rho_freshwater_kg_m3,
                rho_seawater=config.rho_seawater_kg_m3,
            ),
            config.start_month,
            config.end_month,
        ),
    }


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
    wet_global = apply_piecewise_trend_correction(
        raw_global,
        config.wet_tropo_rate_mm_per_year,
        pd.Period(config.wet_tropo_start_month, "M").to_timestamp() + pd.Timedelta(days=14),
        "wet_tropo_corrected",
    )
    corrected_global = apply_scalar_gia(wet_global, config.altimetry_gia_rate_mm_per_year)
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
    wet_common = apply_piecewise_trend_correction(
        raw_common,
        config.wet_tropo_rate_mm_per_year,
        pd.Period(config.wet_tropo_start_month, "M").to_timestamp() + pd.Timedelta(days=14),
        "wet_tropo_corrected",
    )
    corrected_common_base = apply_scalar_gia(wet_common, config.altimetry_gia_rate_mm_per_year)
    corrected_common = MonthlySeries(
        corrected_common_base.time,
        corrected_common_base.values,
        "gmsl_budget_common_wet_gia_corrected",
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
    mass_shc = _renamed(
        mass,
        "ocean_mass_csr_shc_300km",
        role="OBD source consistency diagnostic; not the Jin Mascon ensemble",
    )

    masks: dict[str, SpatialMask] = {
        "altimetry_global": altimetry_mask,
        "budget_common": budget_mask,
    }
    mascon_series: dict[str, MonthlySeries] = {}
    mascon_paths = (config.csr_mascon_path, config.jpl_mascon_path, config.gsfc_sla_mascon_path)
    if any(path is not None for path in mascon_paths):
        spacing = config.mascon_grid_spacing_degrees
        target_latitude = np.arange(-90.0 + spacing / 2.0, 90.0, spacing)
        target_longitude = np.arange(spacing / 2.0, 360.0, spacing)
        mascon_global = _mask_with_hash(
            resample_mask_nearest(altimetry_mask, target_latitude, target_longitude, "mascon_global_1deg"),
            "mascon_global_1deg",
        )
        mascon_300km = _mask_with_hash(
            resample_mask_nearest(budget_mask, target_latitude, target_longitude, "mascon_300km_sensitivity"),
            "mascon_300km_sensitivity",
        )
        masks.update({"mascon_global_1deg": mascon_global, "mascon_300km_sensitivity": mascon_300km})
        observed_global = _read_mascon_centers(config, mascon_global)
        filled_global, ensemble_global = build_ssa_ensemble(
            observed_global, config.start_month, config.end_month
        )
        observed_300km = _read_mascon_centers(config, mascon_300km)
        _, ensemble_300km = build_ssa_ensemble(observed_300km, config.start_month, config.end_month)
        mascon_series.update({item.name: item for item in observed_global.values()})
        mascon_series.update({item.name: item for item in filled_global.values()})
        mascon_series[ensemble_global.name] = ensemble_global
        sensitivity = _renamed(
            ensemble_300km,
            "ocean_mass_ensemble_300km_sensitivity",
            result_role="sensitivity only; not the Jin main domain",
        )
        mascon_series[sensitivity.name] = sensitivity

    input_paths = [
        config.cmems_sla_path,
        config.cdt_land_mask_path,
        config.cdt_distance_path,
        *(item.path for item in epochs),
        *(path for path in mascon_paths if path is not None),
    ]
    provenance = {
        "input_sha256": {str(path.resolve()): sha256_file(path) for path in input_paths},
        "grace_preprocessing_hash": processed.preprocessing_hash,
        "mask_hashes": {
            "altimetry_global": altimetry_mask.metadata["sha256"],
            "budget_common": budget_mask.metadata["sha256"],
            **{name: mask.metadata["sha256"] for name, mask in masks.items() if name.startswith("mascon_")},
        },
        "mascon_gia_policy": "product GIA correction retained; no second GIA correction",
        "mascon_processing_order": "SSA separately for CSR/JPL/GSFC, then arithmetic mean",
    }
    return PipelineRun(
        series={
            item.name: item
            for item in (
                raw_global,
                wet_global,
                corrected_global,
                raw_common,
                wet_common,
                corrected_common,
                mass_shc,
                obd,
            )
        } | mascon_series,
        masks=masks,
        provenance=provenance,
    )


def _monthly_dataset(series: Mapping[str, MonthlySeries]) -> xr.Dataset:
    variables = {}
    for name, item in series.items():
        variables[name] = xr.DataArray(
            item.values,
            dims=("time",),
            coords={"time": ("time", item.time.to_numpy())},
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
    dataset.to_netcdf(path, engine="scipy")


def _coverage_warnings(series: Mapping[str, MonthlySeries], expected_end_month: str) -> list[str]:
    expected = pd.Period(expected_end_month, freq="M")
    warnings = []
    mass_name = "ocean_mass_ensemble" if "ocean_mass_ensemble" in series else "ocean_mass_csr"
    for name in (mass_name, "obd"):
        if name not in series:
            continue
        item = series[name]
        valid = np.isfinite(item.values)
        if not np.any(valid):
            warnings.append(f"{name} has no finite months through requested end month {expected}.")
            continue
        last = item.time[valid].to_period("M").max()
        if last < expected:
            warnings.append(f"{name} ends at {last} before requested end month {expected}.")
    return warnings


def _trend_rows(
    series: Mapping[str, MonthlySeries], hac_lags: int, expected_end_month: str
) -> tuple[list[TrendResult], list[str]]:
    results: list[TrendResult] = []
    warnings = _coverage_warnings(series, expected_end_month)
    for item in series.values():
        try:
            results.append(fit_trend(item, hac_lags))
        except ValueError as error:
            warnings.append(f"Trend unavailable for {item.name}: {error}")
    obd_trend = next((item.trend_mm_per_year for item in results if item.series_name == "obd"), None)
    if obd_trend is not None and not -0.30 <= obd_trend <= 0.10:
        warnings.append(f"OBD trend {obd_trend:.3f} mm/yr is outside the diagnostic range [-0.30, +0.10].")
    return results, warnings


def write_run_outputs(
    run: PipelineRun,
    output_dir: str | Path,
    config: PipelineConfig,
    report_run_id: str | None = None,
) -> None:
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
    dataset.to_netcdf(destination / "monthly_budget.nc", engine="scipy")
    dataset.to_dataframe().to_csv(destination / "monthly_budget.csv", index_label="time")
    trends, warnings = _trend_rows(run.series, config.hac_lags, config.end_month)
    trend_rows = [{**asdict(item), "result_role": "estimated"} for item in trends]
    trend_rows.append(
        {
            "series_name": "gmsl_final_adopted",
            "trend_mm_per_year": config.final_gmsl_trend_mm_per_year,
            "result_role": "adopted_final",
        }
    )
    pd.DataFrame(trend_rows).to_csv(destination / "trend_summary.csv", index=False)
    trend_by_name = {item.series_name: item for item in trends}
    corrected_common = trend_by_name.get("gmsl_budget_common_wet_gia_corrected")
    obd_result = trend_by_name.get("obd")
    computed_budget_trend = (
        corrected_common.trend_mm_per_year + obd_result.trend_mm_per_year
        if corrected_common is not None and obd_result is not None
        else None
    )
    mass_ensemble = trend_by_name.get("ocean_mass_ensemble")
    mass_300km = trend_by_name.get("ocean_mass_ensemble_300km_sensitivity")
    final_summary = {
        "adopted_gmsl_trend_mm_per_year": config.final_gmsl_trend_mm_per_year,
        "adopted_result_role": "final reported value",
        "computed_common_domain_gmsl_plus_obd_mm_per_year": computed_budget_trend,
        "computed_minus_adopted_gmsl_mm_per_year": (
            computed_budget_trend - config.final_gmsl_trend_mm_per_year
            if computed_budget_trend is not None
            else None
        ),
        "ocean_mass_ensemble_trend_mm_per_year": (
            mass_ensemble.trend_mm_per_year if mass_ensemble is not None else None
        ),
        "ocean_mass_ensemble_hac_standard_error_mm_per_year": (
            mass_ensemble.hac_standard_error if mass_ensemble is not None else None
        ),
        "ocean_mass_300km_sensitivity_trend_mm_per_year": (
            mass_300km.trend_mm_per_year if mass_300km is not None else None
        ),
        "jin_ocean_mass_reference_mm_per_year": 1.39,
        "jin_ocean_mass_reference_uncertainty_mm_per_year": 0.32,
        "mascon_domain": "global ocean at 1 degree; no coastal buffer",
        "mascon_300km_role": "sensitivity only",
        "mascon_processing_order": "SSA separately for CSR/JPL/GSFC, then arithmetic mean",
        "mascon_gia_policy": "product GIA correction retained; no second GIA correction",
        "estimated_trends_mm_per_year": {
            name: result.trend_mm_per_year for name, result in trend_by_name.items()
        },
        "note": "The adopted GMSL trend is reported explicitly and is not imposed on any monthly series.",
    }
    write_json(destination / "final_summary.json", final_summary)
    for name, mask in run.masks.items():
        _write_mask(mask, destination / f"mask_{name}.nc")
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
    if "ocean_mass_ensemble" in run.series:
        write_mass_figure(
            run.series,
            destination / "figure1_mass_component.png",
            config.hac_lags,
        )
    write_run_report(
        destination / "run_report.md",
        report_run_id or destination.name,
        config.sha256(),
        trends,
        False,
        warnings,
    )


def run_pipeline(config: PipelineConfig) -> Path:
    config_hash = config.sha256()
    inventory, inventory_hash = _input_inventory(config)
    source_inventory, source_hash = _source_inventory()
    run_id = config.run_id or _default_run_id(config_hash, inventory_hash, source_hash)
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
                and previous.get("processing_source_sha256") == source_hash
            ):
                return destination
            if (
                previous.get("configuration_sha256") == config_hash
                and previous.get("input_inventory_sha256") == inventory_hash
            ):
                raise FileExistsError(f"run directory exists with a different processing source hash: {destination}")
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
            "processing_source_inventory": source_inventory,
            "processing_source_sha256": source_hash,
        },
    )
    staging = output_root / f".{run_id}.tmp-{uuid4().hex}"
    try:
        write_run_outputs(run, staging, config, report_run_id=run_id)
        staging.replace(destination)
    except Exception:
        if staging.exists():
            shutil.rmtree(staging)
        raise
    return destination
