"""Build an auditable four-panel regional-TWS reproduction of Jin Figure 3."""

from __future__ import annotations

import argparse
import csv
from datetime import datetime
import json
from pathlib import Path
from typing import Any, Mapping

import numpy as np

from .figure3.adapters import load_custom_l3, load_mascon
from .figure3.integrate import RegionalSeries, integrate_regions
from .figure3.masks import (
    ContinentMaskSet,
    build_continent_masks,
    read_mask_netcdf,
    write_mask_netcdf,
)
from .figure3.merge import fill_missing_months
from .figure3.metrics import event_metrics
from .figure3.provenance import sha256_file
from .figure3.regrid import nearest_regrid
from .figure3.temporal import (
    ProcessedRegionalSeries,
    RegionalEnsemble,
    combine_centers,
    process_interannual,
)
from .figure3.types import MonthlyGridSeries


_REGION_LABELS = {
    "total": "Total",
    "africa": "Africa",
    "asia": "Asia",
    "europe": "Europe",
    "north_america": "North America",
    "south_america": "South America",
    "oceania": "Oceania",
}
_REGION_COLORS = {
    "total": "#111111",
    "africa": "#d55e00",
    "asia": "#e69f00",
    "europe": "#cc79a7",
    "north_america": "#0072b2",
    "south_america": "#009e73",
    "oceania": "#56b4e9",
}


def _resolve(project_root: Path, value: str | Path) -> Path:
    path = Path(value)
    return path if path.is_absolute() else project_root / path


def _resolved_input_config(config: Mapping[str, Any], project_root: Path) -> dict[str, Any]:
    resolved = dict(config)
    resolved["path"] = str(_resolve(project_root, config["path"]))
    return resolved


def _empty_reconstruction(series: MonthlyGridSeries) -> MonthlyGridSeries:
    return MonthlyGridSeries(
        source_id=f"{series.source_id}-no-reconstruction",
        months=np.asarray([], dtype="U7"),
        lat=series.lat,
        lon=series.lon,
        ewh_mm=np.empty((0, series.lat.size, series.lon.size)),
        valid_month=np.asarray([], dtype=bool),
        month_status=np.asarray([], dtype="U13"),
        metadata={},
    )


def _load_masks(
    config: Mapping[str, Any],
    project_root: Path,
    target_lat: np.ndarray,
    target_lon: np.ndarray,
) -> ContinentMaskSet:
    if "precomputed_path" in config:
        masks = read_mask_netcdf(_resolve(project_root, config["precomputed_path"]))
    else:
        masks = build_continent_masks(
            _resolve(project_root, config["vector_path"]),
            target_lat,
            target_lon,
            coastal_exclusion_km=float(config.get("coastal_exclusion_km", 300.0)),
            coastline_densify_deg=float(config.get("coastline_densify_deg", 0.25)),
            radius_m=float(config.get("earth_radius_m", 6_371_000.0)),
        )
    if not np.array_equal(masks.lat, target_lat) or not np.array_equal(
        masks.lon, target_lon
    ):
        raise ValueError("configured masks do not match the target grid")
    return masks


def _load_centers(
    config: Mapping[str, Any],
    project_root: Path,
    target_lat: np.ndarray,
    target_lon: np.ndarray,
) -> tuple[dict[str, MonthlyGridSeries], dict[str, list[str]], list[Path]]:
    centers: dict[str, MonthlyGridSeries] = {}
    reconstructed_months: dict[str, list[str]] = {}
    input_paths: list[Path] = []
    start = config["time"]["start"]
    end = config["time"]["end"]
    mode = str(config["mode"])

    if mode == "custom_l3":
        center_items = {"custom_l3": {"observed": config["inputs"]["custom_l3"]}}
    else:
        center_items = config["inputs"]["centers"]

    for center_id, center_config in center_items.items():
        observed_config = _resolved_input_config(center_config["observed"], project_root)
        input_paths.append(Path(observed_config["path"]))
        observed = (
            load_custom_l3(observed_config)
            if mode == "custom_l3"
            else load_mascon(observed_config)
        )
        observed = nearest_regrid(observed, target_lat, target_lon)

        reconstruction_config = center_config.get("reconstruction")
        if reconstruction_config is None:
            reconstruction = _empty_reconstruction(observed)
        else:
            resolved_reconstruction = _resolved_input_config(
                reconstruction_config, project_root
            )
            input_paths.append(Path(resolved_reconstruction["path"]))
            reconstruction = nearest_regrid(
                load_mascon(resolved_reconstruction), target_lat, target_lon
            )
        merged = fill_missing_months(
            observed,
            reconstruction,
            start,
            end,
            reconstruction_end=str(
                center_config.get("reconstruction_end", "2022-12")
            ),
        )
        centers[center_id] = merged
        reconstructed_months[center_id] = list(
            merged.metadata.get("reconstructed_months", [])
        )
    return centers, reconstructed_months, input_paths


def _write_plotting_csv(path: Path, ensemble: RegionalEnsemble) -> None:
    fields = ["month"]
    for name in ensemble.region_names:
        fields.extend(
            [
                f"{name}_mean_mm",
                f"{name}_sample_std_mm",
                f"{name}_valid_center_count",
            ]
        )
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for index, month in enumerate(ensemble.months):
            row: dict[str, Any] = {"month": month}
            for name in ensemble.region_names:
                row[f"{name}_mean_mm"] = ensemble.mean_mm[name][index]
                row[f"{name}_sample_std_mm"] = ensemble.sample_std_mm[name][index]
                row[f"{name}_valid_center_count"] = int(
                    ensemble.valid_center_count[name][index]
                )
            writer.writerow(row)


def _write_center_csv(
    path: Path, processed_by_center: Mapping[str, ProcessedRegionalSeries]
) -> None:
    region_names = next(iter(processed_by_center.values())).region_names
    fields = ["center", "month", "month_status"]
    for stage in ("raw", "deseasoned", "detrended", "smoothed"):
        fields.extend(f"{name}_{stage}_mm" for name in region_names)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for center, series in processed_by_center.items():
            for index, month in enumerate(series.months):
                row: dict[str, Any] = {
                    "center": center,
                    "month": month,
                    "month_status": series.month_status[index],
                }
                for name in region_names:
                    row[f"{name}_raw_mm"] = series.raw_mm[name][index]
                    row[f"{name}_deseasoned_mm"] = series.deseasoned_mm[name][index]
                    row[f"{name}_detrended_mm"] = series.detrended_mm[name][index]
                    row[f"{name}_smoothed_mm"] = series.smoothed_mm[name][index]
                writer.writerow(row)


def _write_metrics_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]) if rows else [])
        if rows:
            writer.writeheader()
            writer.writerows(rows)


def _plot(
    ensemble: RegionalEnsemble,
    metrics_rows: list[dict[str, Any]],
    events: list[Mapping[str, str]],
    png: Path,
    pdf: Path,
) -> None:
    import matplotlib.pyplot as plt
    import matplotlib.dates as mdates

    if len(events) != 2:
        raise ValueError("Figure 3 requires exactly two configured events")
    dates = np.asarray(
        [datetime(int(month[:4]), int(month[5:7]), 15) for month in ensemble.months]
    )
    fig, axes = plt.subplots(2, 2, figsize=(12.4, 7.6))
    panel_letters = (("a", "b"), ("c", "d"))
    order = ("total",) + tuple(
        name for name in ensemble.region_names if name != "total"
    )
    for row_index, event in enumerate(events):
        selection = (ensemble.months >= event["display_start"]) & (
            ensemble.months <= event["display_end"]
        )
        line_axis = axes[row_index, 0]
        for name in order:
            line_axis.plot(
                dates[selection],
                ensemble.mean_mm[name][selection],
                label=_REGION_LABELS[name],
                color=_REGION_COLORS[name],
                linewidth=2.2 if name == "total" else 1.25,
            )
        line_axis.axhline(0.0, color="#777777", linewidth=0.7)
        line_axis.axvspan(
            datetime(int(event["start"][:4]), int(event["start"][5:7]), 1),
            datetime(int(event["end"][:4]), int(event["end"][5:7]), 28),
            color="#dddddd",
            alpha=0.55,
            zorder=0,
        )
        line_axis.set_title(event.get("title", event["id"]))
        line_axis.set_ylabel("Regional contribution to GMSL (mm)")
        line_axis.grid(True, alpha=0.2)
        line_axis.xaxis.set_major_locator(mdates.AutoDateLocator(minticks=4, maxticks=7))
        line_axis.xaxis.set_major_formatter(mdates.DateFormatter("%Y-%m"))
        line_axis.text(
            0.02,
            0.95,
            f"({panel_letters[row_index][0]})",
            transform=line_axis.transAxes,
            va="top",
            fontweight="bold",
        )
        if row_index == 0:
            line_axis.legend(ncol=2, fontsize=7.3, frameon=False)

        bar_axis = axes[row_index, 1]
        event_rows = {row["region"]: row for row in metrics_rows if row["event_id"] == event["id"]}
        changes = [event_rows[name]["change_mm"] for name in order]
        x = np.arange(len(order))
        bar_axis.bar(x, changes, color=[_REGION_COLORS[name] for name in order])
        bar_axis.axhline(0.0, color="#777777", linewidth=0.7)
        bar_axis.set_xticks(x)
        bar_axis.set_xticklabels([_REGION_LABELS[name] for name in order], rotation=24, ha="right")
        bar_axis.set_ylabel("Endpoint change (mm)")
        bar_axis.grid(True, axis="y", alpha=0.2)
        bar_axis.text(
            0.02,
            0.95,
            f"({panel_letters[row_index][1]})",
            transform=bar_axis.transAxes,
            va="top",
            fontweight="bold",
        )
    fig.suptitle("Regional terrestrial-water-storage contributions during two El Niño events")
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    fig.savefig(png, dpi=300)
    fig.savefig(pdf)
    plt.close(fig)


def _event_window_gaps(
    ensemble: RegionalEnsemble, events: list[Mapping[str, str]]
) -> list[dict[str, str]]:
    gaps: list[dict[str, str]] = []
    for event in events:
        selection = (ensemble.months >= event["start"]) & (
            ensemble.months <= event["end"]
        )
        for name in ensemble.region_names:
            for month in ensemble.months[selection][
                ~np.isfinite(ensemble.mean_mm[name][selection])
            ]:
                gaps.append({"event_id": event["id"], "region": name, "month": month})
    return gaps


def build_figure3(
    config_path: Path | str,
    project_root: Path | str | None = None,
) -> dict[str, Path]:
    """Run the configured Figure 3 workflow and return all artifact paths."""
    config_path = Path(config_path)
    root = (
        Path(project_root)
        if project_root is not None
        else Path(__file__).resolve().parents[2]
    )
    if not config_path.is_absolute():
        config_path = root / config_path
    config = json.loads(config_path.read_text(encoding="utf-8"))
    corrections = config.get("corrections", {})
    if corrections.get("apply_gia") or corrections.get("apply_obd"):
        raise ValueError("Figure 3 regional TWS pipeline applies neither extra GIA nor OBD")

    target_lat = np.asarray(config["target_grid"]["lat"], dtype=float)
    target_lon = np.asarray(config["target_grid"]["lon"], dtype=float)
    masks = _load_masks(config["mask"], root, target_lat, target_lon)
    centers, reconstructed_months, input_paths = _load_centers(
        config, root, target_lat, target_lon
    )
    regional_by_center: dict[str, RegionalSeries] = {
        center: integrate_regions(
            series,
            masks,
            ocean_area_m2=float(config["integration"]["global_ocean_area_m2"]),
        )
        for center, series in centers.items()
    }
    processed = {
        center: process_interannual(series)
        for center, series in regional_by_center.items()
    }
    events = list(config["events"])
    paper_mode = str(config["mode"]).startswith("paper")
    ensemble = combine_centers(
        processed, paper_mode=paper_mode, events=events
    )
    metrics_rows = event_metrics(
        ensemble, events, config.get("paper_references", {})
    )

    output_directory = _resolve(root, config["output"]["directory"])
    output_directory.mkdir(parents=True, exist_ok=True)
    stem = config["output"].get("stem", "figure03_regional_tws")
    outputs = {
        "png": output_directory / f"{stem}.png",
        "pdf": output_directory / f"{stem}.pdf",
        "plotting_data": output_directory / f"{stem}_plotting_data.csv",
        "regional_by_center": output_directory / f"{stem}_regional_by_center.csv",
        "metrics": output_directory / f"{stem}_metrics.csv",
        "masks_netcdf": output_directory / f"{stem}_masks.nc",
        "config_snapshot": output_directory / f"{stem}_config.json",
        "method_report": output_directory / f"{stem}_method.md",
        "manifest": output_directory / f"{stem}_manifest.json",
    }
    _write_plotting_csv(outputs["plotting_data"], ensemble)
    _write_center_csv(outputs["regional_by_center"], processed)
    _write_metrics_csv(outputs["metrics"], metrics_rows)
    write_mask_netcdf(masks, outputs["masks_netcdf"])
    outputs["config_snapshot"].write_text(
        json.dumps(config, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    outputs["method_report"].write_text(
        "# Figure 3 regional TWS method\n\n"
        "Each center is normalized and regridded separately. Original observed "
        "months take precedence; registered Xie–Yi values fill missing months only. "
        "Six Natural Earth continental domains exclude Greenland, Antarctica, and "
        "land within 300 km of the exterior ocean coastline. Regional EWH mass is "
        "divided by the configured global ocean area. Each center is deseasoned, "
        "OLS-detrended, and strictly 3-month centered-smoothed before the arithmetic "
        "center mean is formed. No OBD and no additional GIA correction are applied.\n",
        encoding="utf-8",
    )
    _plot(ensemble, metrics_rows, events, outputs["png"], outputs["pdf"])

    manifest = {
        "mode": config["mode"],
        "config_path": str(config_path),
        "config_sha256": sha256_file(config_path),
        "input_sha256": {
            str(path): sha256_file(path) for path in sorted(set(input_paths))
        },
        "reconstructed_months_by_center": {
            center: reconstructed_months[center] for center in sorted(reconstructed_months)
        },
        "event_window_gaps": _event_window_gaps(ensemble, events),
        "corrections": {
            "gia_applied_in_pipeline": False,
            "obd_applied": False,
        },
        "ocean_area_m2": float(config["integration"]["global_ocean_area_m2"]),
        "processing_order": [
            "monthly_climatology",
            "ols_detrend",
            "centered_3_month_mean",
            "center_arithmetic_mean",
        ],
        "center_ids": list(ensemble.center_ids),
        "metrics": metrics_rows,
        "output_sha256": {
            key: sha256_file(path)
            for key, path in outputs.items()
            if key != "manifest"
        },
    }
    outputs["manifest"].write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return outputs


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True, type=Path)
    arguments = parser.parse_args()
    outputs = build_figure3(arguments.config)
    for name, path in outputs.items():
        print(f"{name}: {path}")


if __name__ == "__main__":
    main()
