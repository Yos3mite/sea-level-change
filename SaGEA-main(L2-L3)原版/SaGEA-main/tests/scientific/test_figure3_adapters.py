from pathlib import Path
import warnings

import numpy as np
import pytest
from netCDF4 import Dataset

from pysrc.reference_products.figure3.adapters import load_custom_l3, load_mascon
from pysrc.reference_products.figure3.types import MonthlyGridSeries


def _write_grid_file(
    path: Path,
    *,
    units: str = "cm",
    include_valid_month: bool = False,
) -> Path:
    with Dataset(path, "w") as dataset:
        dataset.createDimension("time", 2)
        dataset.createDimension("lat", 2)
        dataset.createDimension("lon", 3)

        time = dataset.createVariable("time", "f8", ("time",))
        time.units = "days since 2023-01-01 00:00:00"
        time.calendar = "standard"
        time[:] = [14.0, 45.0]

        lat = dataset.createVariable("latitude", "f8", ("lat",))
        lat[:] = [0.5, -0.5]
        lon = dataset.createVariable("longitude", "f8", ("lon",))
        lon[:] = [359.5, 0.5, 1.5]

        field = dataset.createVariable(
            "lwe",
            "f8",
            ("lon", "lat", "time"),
            fill_value=-9999.0,
        )
        field.units = units
        values = np.full((3, 2, 2), 2.0)
        values[0, 1, 0] = 1.0
        values[2, 0, 1] = -9999.0
        with warnings.catch_warnings():
            warnings.filterwarnings(
                "ignore",
                message="Setting the shape on a NumPy array has been deprecated",
                category=DeprecationWarning,
            )
            field[:, :, :] = values

        if include_valid_month:
            valid = dataset.createVariable("valid_month", "i1", ("time",))
            valid[:] = [1, 0]
    return path


def _config(path: Path, *, status: str = "observed") -> dict:
    return {
        "source_id": "test-grid",
        "path": str(path),
        "variables": {
            "time": "time",
            "lat": "latitude",
            "lon": "longitude",
            "field": "lwe",
        },
        "dimensions": {"time": "time", "lat": "lat", "lon": "lon"},
        "month_status": status,
        "metadata": {"gia_corrected": True},
    }


def test_monthly_grid_series_rejects_duplicate_months():
    """Catch products that could silently overwrite one of two monthly fields."""
    with pytest.raises(ValueError, match="duplicate months"):
        MonthlyGridSeries(
            source_id="x",
            months=np.array(["2023-01", "2023-01"]),
            lat=np.array([-0.5, 0.5]),
            lon=np.array([-0.5, 0.5]),
            ewh_mm=np.zeros((2, 2, 2)),
            valid_month=np.ones(2, dtype=bool),
            month_status=np.array(["observed", "observed"]),
            metadata={},
        )


def test_monthly_grid_series_requires_missing_month_fields_to_be_nan():
    """Catch contradictory records marked missing while still carrying data."""
    with pytest.raises(ValueError, match="missing months must contain only NaN"):
        MonthlyGridSeries(
            source_id="x",
            months=np.array(["2023-01"]),
            lat=np.array([-0.5, 0.5]),
            lon=np.array([-0.5, 0.5]),
            ewh_mm=np.zeros((1, 2, 2)),
            valid_month=np.array([False]),
            month_status=np.array(["missing"]),
            metadata={},
        )


def test_load_mascon_normalizes_units_axis_order_and_coordinates(tmp_path: Path):
    """Catch cm fields or descending/wrapped coordinates leaking downstream."""
    path = _write_grid_file(tmp_path / "mascon.nc")

    series = load_mascon(_config(path))

    assert series.ewh_mm.shape == (2, 2, 3)
    assert series.ewh_mm[0, 0, 0] == pytest.approx(10.0)
    assert np.isnan(series.ewh_mm[1, 1, 2])
    assert series.months.tolist() == ["2023-01", "2023-02"]
    assert series.lat.tolist() == [-0.5, 0.5]
    assert series.lon.tolist() == [-0.5, 0.5, 1.5]
    assert series.valid_month.tolist() == [True, True]
    assert series.month_status.tolist() == ["observed", "observed"]
    assert series.metadata["gia_corrected"] is True
    assert series.metadata["source_units"] == "cm"


def test_load_custom_l3_respects_valid_month_and_masks_invalid_field(tmp_path: Path):
    """Catch local Level-3 validity flags being ignored during normalization."""
    path = _write_grid_file(
        tmp_path / "custom.nc",
        units="m",
        include_valid_month=True,
    )
    config = _config(path)
    config["variables"]["valid_month"] = "valid_month"

    series = load_custom_l3(config)

    assert series.valid_month.tolist() == [True, False]
    assert series.month_status.tolist() == ["observed", "missing"]
    assert series.ewh_mm[0, 0, 0] == pytest.approx(1000.0)
    assert np.isnan(series.ewh_mm[1]).all()


def test_load_mascon_rejects_unrecognized_water_height_units(tmp_path: Path):
    """Catch accidental ingestion of mass, pressure, or unitless fields as EWH."""
    path = _write_grid_file(tmp_path / "unknown.nc", units="kg m-2")

    with pytest.raises(ValueError, match="unsupported water-height units"):
        load_mascon(_config(path, status="reconstructed"))


def test_load_mascon_can_select_a_requested_month_range(tmp_path: Path):
    """Catch full multi-decade grids being loaded when only the study window is needed."""
    path = _write_grid_file(tmp_path / "windowed.nc")
    config = _config(path)
    config["time_start"] = "2023-02"
    config["time_end"] = "2023-02"

    series = load_mascon(config)

    assert series.months.tolist() == ["2023-02"]
    assert series.ewh_mm.shape == (1, 2, 3)


def test_load_custom_l3_accepts_string_months_and_mm_ewh_units(tmp_path: Path):
    """Catch the project's native custom Level-3 time/unit contract being rejected."""
    path = tmp_path / "native_custom.nc"
    with Dataset(path, "w") as dataset:
        dataset.createDimension("time", 2)
        dataset.createDimension("latitude", 1)
        dataset.createDimension("longitude", 1)
        time = dataset.createVariable("time", str, ("time",))
        time[:] = np.asarray(["2023-01", "2023-02"], dtype=object)
        dataset.createVariable("lat", "f8", ("latitude",))[:] = [0.0]
        dataset.createVariable("lon", "f8", ("longitude",))[:] = [0.0]
        field = dataset.createVariable("field", "f8", ("time", "latitude", "longitude"))
        field.units = "mm EWH"
        with warnings.catch_warnings():
            warnings.filterwarnings(
                "ignore",
                message="Setting the shape on a NumPy array has been deprecated",
                category=DeprecationWarning,
            )
            field[:] = np.asarray([[[1.0]], [[2.0]]])
        dataset.createVariable("valid_month", "i1", ("time",))[:] = [1, 1]
    config = {
        "source_id": "native-custom",
        "path": str(path),
        "variables": {
            "time": "time",
            "lat": "lat",
            "lon": "lon",
            "field": "field",
            "valid_month": "valid_month",
        },
        "dimensions": {
            "time": "time",
            "lat": "latitude",
            "lon": "longitude",
        },
    }

    series = load_custom_l3(config)

    assert series.months.tolist() == ["2023-01", "2023-02"]
    assert series.ewh_mm[:, 0, 0].tolist() == [1.0, 2.0]
