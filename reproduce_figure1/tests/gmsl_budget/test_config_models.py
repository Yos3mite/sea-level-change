import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from gmsl_budget.config import PipelineConfig
from gmsl_budget.models import MonthlySeries, SpatialMask


def _config_payload(**overrides):
    payload = {
        "cmems_sla_path": "D:/temp_sealevel_data/cmems_sla_monthly.nc",
        "cmems_indicator_path": "D:/temp_sealevel_data/cmems_gmsl_indicator.nc",
        "cdt_land_mask_path": "D:/data/land_mask.mat",
        "cdt_distance_path": "D:/data/distance2coast.mat",
        "sagea_root": "D:/AAAA海平面变化/SaGEA-main(L2-L3)/SaGEA-main",
        "grace_gfc_dir": "D:/temp_sealevel_data/grace_csr_ddk1",
        "output_root": "D:/AAAA海平面变化/reproduce_figure1/output/optimized_budget",
        "start_month": "2013-11",
        "end_month": "2024-10",
        "gia_mode": "scalar",
        "altimetry_gia_rate_mm_per_year": 0.30,
        "rho_freshwater_kg_m3": 1000.0,
        "rho_seawater_kg_m3": 1028.0,
        "min_valid_weight_fraction": 0.995,
        "hac_lags": 12,
        "coastal_buffer_km": 300.0,
        "grid_spacing_degrees": 0.25,
        "mass_preprocessing": {"release": "RL06.3", "filter": "DDK1", "gia": "Caron2018"},
        "obd_preprocessing": {"release": "RL06.3", "filter": "DDK1", "gia": "Caron2018"},
        "steric_path": None,
        "run_id": None,
    }
    payload.update(overrides)
    return payload


def _write_config(directory: Path, payload: dict) -> Path:
    directory.mkdir(parents=True, exist_ok=True)
    path = directory / "config.json"
    path.write_text(json.dumps(payload), encoding="utf-8")
    return path


def test_config_rejects_non_positive_density(tmp_path):
    path = _write_config(tmp_path, _config_payload(rho_seawater_kg_m3=0.0))
    with pytest.raises(ValueError, match="rho_seawater_kg_m3"):
        PipelineConfig.load(path)


def test_config_rejects_different_mass_and_obd_preprocessing(tmp_path):
    path = _write_config(
        tmp_path,
        _config_payload(obd_preprocessing={"release": "RL06.3", "filter": "Gaussian", "gia": "Caron2018"}),
    )
    with pytest.raises(ValueError, match="identical preprocessing"):
        PipelineConfig.load(path)


def test_config_hash_is_independent_of_json_key_order(tmp_path):
    payload = _config_payload()
    left = PipelineConfig.load(_write_config(tmp_path / "left", payload))
    reversed_payload = dict(reversed(list(payload.items())))
    right = PipelineConfig.load(_write_config(tmp_path / "right", reversed_payload))
    assert left.sha256() == right.sha256()


def test_monthly_series_rejects_duplicate_months():
    time = pd.to_datetime(["2020-01-15", "2020-01-20"])
    with pytest.raises(ValueError, match="duplicate month"):
        MonthlySeries(time, np.array([1.0, 2.0]), "x", "mm", {})


def test_monthly_series_normalizes_times_to_month_midpoint():
    result = MonthlySeries(
        pd.to_datetime(["2020-01-01", "2020-02-29"]),
        np.array([1.0, 2.0]),
        "x",
        "mm",
        {},
    )
    assert result.time.tolist() == pd.to_datetime(["2020-01-15", "2020-02-15"]).tolist()
    assert result.values.dtype == np.float64


def test_spatial_mask_rejects_fraction_outside_zero_one():
    with pytest.raises(ValueError, match="ocean_fraction"):
        SpatialMask(
            latitude=np.array([0.0]),
            longitude=np.array([0.0, 1.0]),
            ocean_fraction=np.array([[0.5, 1.1]]),
            support=np.array([[True, True]]),
            metadata={},
        )
