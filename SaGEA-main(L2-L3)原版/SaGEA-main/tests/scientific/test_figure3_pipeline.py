from datetime import datetime
import json
from pathlib import Path
import warnings

from netCDF4 import Dataset, date2num
import numpy as np
import pytest

from pysrc.reference_products.build_figure3_regional_tws import build_figure3
from pysrc.reference_products.figure3.masks import ContinentMaskSet, write_mask_netcdf


pytestmark = pytest.mark.filterwarnings("ignore::DeprecationWarning:matplotlib\\..*")


def _month_axis(start: str, count: int) -> list[str]:
    first = np.datetime64(start, "M")
    return [str(item) for item in first + np.arange(count)]


def _write_center(
    path: Path,
    months: list[str],
    *,
    center_offset: float,
    missing_month: str | None = None,
) -> Path:
    with Dataset(path, "w") as dataset:
        dataset.createDimension("time", len(months))
        dataset.createDimension("lat", 1)
        dataset.createDimension("lon", 6)
        time = dataset.createVariable("time", "f8", ("time",))
        time.units = "days since 2019-01-01 00:00:00"
        dates = [datetime(int(month[:4]), int(month[5:7]), 15) for month in months]
        time[:] = date2num(dates, time.units)
        dataset.createVariable("lat", "f8", ("lat",))[:] = [0.0]
        dataset.createVariable("lon", "f8", ("lon",))[:] = np.arange(6)
        field = dataset.createVariable(
            "ewh", "f8", ("time", "lat", "lon"), fill_value=-9999.0
        )
        field.units = "mm"
        valid = dataset.createVariable("valid_month", "i1", ("time",))
        values = np.empty((len(months), 1, 6), dtype=float)
        flags = np.ones(len(months), dtype=np.int8)
        for index, month in enumerate(months):
            values[index, 0] = center_offset + index + np.arange(6)
            if month == missing_month:
                values[index] = -9999.0
                flags[index] = 0
        with warnings.catch_warnings():
            warnings.filterwarnings(
                "ignore",
                message="Setting the shape on a NumPy array has been deprecated",
                category=DeprecationWarning,
            )
            field[:] = values
        valid[:] = flags
    return path


def _input_config(path: Path, source_id: str, status: str) -> dict:
    return {
        "source_id": source_id,
        "path": str(path),
        "variables": {
            "time": "time",
            "lat": "lat",
            "lon": "lon",
            "field": "ewh",
            "valid_month": "valid_month",
        },
        "dimensions": {"time": "time", "lat": "lat", "lon": "lon"},
        "month_status": status,
        "metadata": {"gia_corrected": True},
    }


def test_build_figure3_creates_auditable_artifact_bundle(tmp_path: Path):
    """Catch orchestration that omits reconstructed-month or correction provenance."""
    months = _month_axis("2019-01", 24)
    reconstruction_month = "2020-06"
    centers: dict[str, dict] = {}
    for center_index, center in enumerate(("CSR", "JPL", "GSFC")):
        observed = _write_center(
            tmp_path / f"{center}_observed.nc",
            months,
            center_offset=float(center_index),
            missing_month=reconstruction_month,
        )
        reconstruction = _write_center(
            tmp_path / f"{center}_reconstruction.nc",
            [reconstruction_month],
            center_offset=float(center_index + 17),
        )
        centers[center] = {
            "observed": _input_config(observed, center, "observed"),
            "reconstruction": _input_config(
                reconstruction, f"{center}-Xie-Yi", "reconstructed"
            ),
        }

    masks = ContinentMaskSet(
        lat=np.asarray([0.0]),
        lon=np.arange(6, dtype=float),
        region_names=(
            "africa",
            "asia",
            "europe",
            "north_america",
            "south_america",
            "oceania",
        ),
        region_id=np.asarray([[1, 2, 3, 4, 5, 6]], dtype=np.int16),
        land_mask=np.ones((1, 6), dtype=bool),
        distance_to_coast_km=np.full((1, 6), 500.0),
        coastal_buffer_excluded=np.zeros((1, 6), dtype=bool),
        cell_area_m2=np.ones((1, 6), dtype=float),
        metadata={"fixture": True},
    )
    mask_path = write_mask_netcdf(masks, tmp_path / "mask.nc")
    output_directory = tmp_path / "outputs"
    config = {
        "mode": "paper_mascon",
        "inputs": {"centers": centers},
        "target_grid": {"lat": [0.0], "lon": [0, 1, 2, 3, 4, 5]},
        "mask": {"precomputed_path": str(mask_path)},
        "integration": {"global_ocean_area_m2": 6.0},
        "time": {"start": "2019-01", "end": "2020-12"},
        "events": [
            {
                "id": "event_2019",
                "title": "Event 2019",
                "display_start": "2019-02",
                "display_end": "2019-11",
                "start": "2019-03",
                "end": "2019-10",
            },
            {
                "id": "event_2020",
                "title": "Event 2020",
                "display_start": "2020-02",
                "display_end": "2020-11",
                "start": "2020-03",
                "end": "2020-10",
            },
        ],
        "paper_references": {"event_2019": {}, "event_2020": {}},
        "corrections": {"apply_gia": False, "apply_obd": False},
        "output": {"directory": str(output_directory), "stem": "figure03_test"},
    }
    config_path = tmp_path / "figure3.json"
    config_path.write_text(json.dumps(config), encoding="utf-8")

    outputs = build_figure3(config_path, project_root=tmp_path)

    expected = {
        "png",
        "pdf",
        "plotting_data",
        "regional_by_center",
        "metrics",
        "masks_netcdf",
        "config_snapshot",
        "method_report",
        "manifest",
    }
    assert set(outputs) == expected
    assert all(path.is_file() for path in outputs.values())
    manifest = json.loads(outputs["manifest"].read_text(encoding="utf-8"))
    assert manifest["reconstructed_months_by_center"] == {
        "CSR": [reconstruction_month],
        "GSFC": [reconstruction_month],
        "JPL": [reconstruction_month],
    }
    assert manifest["event_window_gaps"] == []
    assert manifest["corrections"] == {
        "gia_applied_in_pipeline": False,
        "obd_applied": False,
    }
