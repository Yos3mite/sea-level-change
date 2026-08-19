import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
import xarray as xr

from gmsl_budget.config import PipelineConfig
from gmsl_budget.models import MonthlySeries, SpatialMask
from gmsl_budget.pipeline import PipelineRun, run_pipeline
from gmsl_budget.provenance import sha256_file


def _config(tmp_path: Path, gia_rate=0.30, run_id="fixture-run") -> PipelineConfig:
    payload = {
        "cmems_sla_path": str(tmp_path / "sla.nc"),
        "cmems_indicator_path": None,
        "cdt_land_mask_path": str(tmp_path / "land.mat"),
        "cdt_distance_path": str(tmp_path / "distance.mat"),
        "sagea_root": str(tmp_path / "sagea"),
        "grace_gfc_dir": str(tmp_path / "gfc"),
        "output_root": str(tmp_path / "output"),
        "start_month": "2018-01",
        "end_month": "2021-12",
        "gia_mode": "scalar",
        "altimetry_gia_rate_mm_per_year": gia_rate,
        "rho_freshwater_kg_m3": 1000.0,
        "rho_seawater_kg_m3": 1028.0,
        "min_valid_weight_fraction": 0.995,
        "hac_lags": 12,
        "coastal_buffer_km": 300.0,
        "grid_spacing_degrees": 0.25,
        "mass_preprocessing": {"filter": "DDK1"},
        "obd_preprocessing": {"filter": "DDK1"},
        "steric_path": None,
        "run_id": run_id,
    }
    path = tmp_path / f"config-{gia_rate}.json"
    path.write_text(json.dumps(payload), encoding="utf-8")
    return PipelineConfig.load(path)


def _fixture_run(config: PipelineConfig) -> PipelineRun:
    time = pd.date_range("2018-01-15", periods=48, freq="MS") + pd.Timedelta(days=14)
    year = np.arange(48, dtype=np.float64) / 12.0
    metadata = {"mask_hash": "budget-mask", "preprocessing_hash": "shared"}
    series = {
        "gmsl_raw": MonthlySeries(time, 4.0 * year, "gmsl_raw", "mm", {"gia_corrected": False}),
        "gmsl_gia_corrected": MonthlySeries(
            time, 4.3 * year, "gmsl_gia_corrected", "mm", {"gia_corrected": True}
        ),
        "ocean_mass_csr": MonthlySeries(time, 2.0 * year, "ocean_mass_csr", "mm", metadata),
        "obd": MonthlySeries(time, -0.1 * year, "obd", "mm", metadata),
    }
    mask = SpatialMask(
        np.array([-45.0, 45.0]),
        np.array([0.0, 180.0]),
        np.ones((2, 2)),
        np.ones((2, 2), dtype=bool),
        {"sha256": "budget-mask", "name": "budget_common"},
    )
    return PipelineRun(series=series, masks={"budget_common": mask, "altimetry_global": mask}, provenance={})


def test_sha256_file_uses_file_bytes(tmp_path):
    path = tmp_path / "input.bin"
    path.write_bytes(b"sea-level")
    assert sha256_file(path) == "f53c4fccbe7d1a73d53fb1551e76725b7319cf7097d78b72a4bee9a9ca365471"


def test_pipeline_without_steric_writes_valid_partial_budget(tmp_path, monkeypatch):
    config = _config(tmp_path)
    monkeypatch.setattr("gmsl_budget.pipeline._compute_run", _fixture_run)
    run_dir = run_pipeline(config)
    with xr.open_dataset(run_dir / "monthly_budget.nc") as dataset:
        assert {"gmsl_raw", "gmsl_gia_corrected", "ocean_mass_csr", "obd"} <= set(dataset)
        assert "steric" not in dataset
        assert "closure" not in dataset
    provenance = json.loads((run_dir / "provenance.json").read_text(encoding="utf-8"))
    assert provenance["budget_closure_available"] is False
    assert (run_dir / "diagnostics.png").is_file()
    assert (run_dir / "run_report.md").is_file()


def test_pipeline_refuses_to_overwrite_different_config(tmp_path, monkeypatch):
    first = _config(tmp_path, gia_rate=0.30)
    second = _config(tmp_path, gia_rate=0.25)
    monkeypatch.setattr("gmsl_budget.pipeline._compute_run", _fixture_run)
    run_pipeline(first)
    with pytest.raises(FileExistsError, match="different configuration hash"):
        run_pipeline(second)
