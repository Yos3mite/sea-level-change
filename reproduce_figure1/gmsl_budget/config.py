from __future__ import annotations

from dataclasses import asdict, dataclass
from hashlib import sha256
import json
import math
from pathlib import Path
from typing import Any, Mapping

import pandas as pd


@dataclass(frozen=True)
class PipelineConfig:
    cmems_sla_path: Path
    cmems_indicator_path: Path | None
    cdt_land_mask_path: Path
    cdt_distance_path: Path
    sagea_root: Path
    grace_gfc_dir: Path
    output_root: Path
    start_month: str
    end_month: str
    gia_mode: str
    altimetry_gia_rate_mm_per_year: float
    rho_freshwater_kg_m3: float
    rho_seawater_kg_m3: float
    min_valid_weight_fraction: float
    hac_lags: int
    coastal_buffer_km: float
    grid_spacing_degrees: float
    mass_preprocessing: Mapping[str, Any]
    obd_preprocessing: Mapping[str, Any]
    steric_path: Path | None = None
    run_id: str | None = None
    csr_mascon_path: Path | None = None
    jpl_mascon_path: Path | None = None
    gsfc_sla_mascon_path: Path | None = None
    wet_tropo_rate_mm_per_year: float = -0.50
    wet_tropo_start_month: str = "2016-01"
    final_gmsl_trend_mm_per_year: float = 3.057
    mascon_grid_spacing_degrees: float = 1.0

    @classmethod
    def load(cls, path: str | Path) -> "PipelineConfig":
        config_path = Path(path)
        with config_path.open("r", encoding="utf-8") as stream:
            payload = json.load(stream)
        path_fields = (
            "cmems_sla_path",
            "cdt_land_mask_path",
            "cdt_distance_path",
            "sagea_root",
            "grace_gfc_dir",
            "output_root",
        )
        for field in path_fields:
            if field not in payload:
                raise ValueError(f"missing required configuration field: {field}")
            payload[field] = Path(payload[field])
        for field in (
            "cmems_indicator_path",
            "steric_path",
            "csr_mascon_path",
            "jpl_mascon_path",
            "gsfc_sla_mascon_path",
        ):
            if payload.get(field) is not None:
                payload[field] = Path(payload[field])
        config = cls(**payload)
        config._validate()
        return config

    def _validate(self) -> None:
        if self.rho_freshwater_kg_m3 <= 0:
            raise ValueError("rho_freshwater_kg_m3 must be positive")
        if self.rho_seawater_kg_m3 <= 0:
            raise ValueError("rho_seawater_kg_m3 must be positive")
        if not 0 < self.min_valid_weight_fraction <= 1:
            raise ValueError("min_valid_weight_fraction must be in (0, 1]")
        if self.gia_mode not in {"scalar", "spatial_caron"}:
            raise ValueError("gia_mode must be scalar or spatial_caron")
        if self.hac_lags < 0:
            raise ValueError("hac_lags must be non-negative")
        if self.coastal_buffer_km < 0:
            raise ValueError("coastal_buffer_km must be non-negative")
        if self.grid_spacing_degrees <= 0:
            raise ValueError("grid_spacing_degrees must be positive")
        if self.mascon_grid_spacing_degrees <= 0:
            raise ValueError("mascon_grid_spacing_degrees must be positive")
        if not math.isfinite(self.altimetry_gia_rate_mm_per_year):
            raise ValueError("altimetry_gia_rate_mm_per_year must be finite")
        if not math.isfinite(self.wet_tropo_rate_mm_per_year):
            raise ValueError("wet_tropo_rate_mm_per_year must be finite")
        if not math.isfinite(self.final_gmsl_trend_mm_per_year):
            raise ValueError("final_gmsl_trend_mm_per_year must be finite")
        start = pd.Period(self.start_month, freq="M")
        end = pd.Period(self.end_month, freq="M")
        wet_tropo_start = pd.Period(self.wet_tropo_start_month, freq="M")
        if start > end:
            raise ValueError("start_month must not be after end_month")
        if wet_tropo_start > end:
            raise ValueError("wet_tropo_start_month must not be after end_month")
        if dict(self.mass_preprocessing) != dict(self.obd_preprocessing):
            raise ValueError("mass and OBD must use identical preprocessing")

    def resolved_dict(self) -> dict[str, Any]:
        result = asdict(self)
        for key, value in tuple(result.items()):
            if isinstance(value, Path):
                result[key] = str(value.resolve())
        return result

    def sha256(self) -> str:
        serialized = json.dumps(
            self.resolved_dict(),
            ensure_ascii=False,
            sort_keys=True,
            separators=(",", ":"),
        ).encode("utf-8")
        return sha256(serialized).hexdigest()
