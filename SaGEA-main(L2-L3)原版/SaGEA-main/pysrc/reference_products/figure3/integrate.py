"""Convert gridded terrestrial EWH into regional global-ESL contributions."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Mapping

import numpy as np

from .masks import ContinentMaskSet
from .types import MonthlyGridSeries


@dataclass(frozen=True)
class RegionalSeries:
    """Per-region monthly contributions expressed as global mean sea-level mm."""

    source_id: str
    months: np.ndarray
    region_names: tuple[str, ...]
    values_mm: Mapping[str, np.ndarray]
    contributing_area_m2: Mapping[str, np.ndarray]
    valid_cell_fraction: Mapping[str, np.ndarray]
    month_status: np.ndarray
    metadata: Mapping[str, Any]

    def __post_init__(self) -> None:
        months = np.asarray(self.months, dtype="U7")
        statuses = np.asarray(self.month_status, dtype="U13")
        if statuses.shape != (months.size,):
            raise ValueError("month_status must contain one value per month")
        expected_names = tuple(self.region_names)
        for collection_name in (
            "values_mm",
            "contributing_area_m2",
            "valid_cell_fraction",
        ):
            collection = getattr(self, collection_name)
            if set(collection) != set(expected_names):
                raise ValueError(f"{collection_name} must contain every region")
            if any(np.asarray(values).shape != (months.size,) for values in collection.values()):
                raise ValueError(f"{collection_name} arrays must match the month axis")
        object.__setattr__(self, "months", months)
        object.__setattr__(self, "region_names", expected_names)
        object.__setattr__(self, "month_status", statuses)
        object.__setattr__(
            self,
            "values_mm",
            {name: np.asarray(value, dtype=float) for name, value in self.values_mm.items()},
        )
        object.__setattr__(
            self,
            "contributing_area_m2",
            {
                name: np.asarray(value, dtype=float)
                for name, value in self.contributing_area_m2.items()
            },
        )
        object.__setattr__(
            self,
            "valid_cell_fraction",
            {
                name: np.asarray(value, dtype=float)
                for name, value in self.valid_cell_fraction.items()
            },
        )
        object.__setattr__(self, "metadata", dict(self.metadata))


def integrate_regions(
    series: MonthlyGridSeries,
    masks: ContinentMaskSet,
    ocean_area_m2: float,
    *,
    water_density_kg_m3: float = 1000.0,
) -> RegionalSeries:
    """Area-integrate regional EWH and normalize by one global ocean area."""
    if ocean_area_m2 <= 0 or water_density_kg_m3 <= 0:
        raise ValueError("ocean area and water density must be positive")
    if not np.array_equal(series.lat, masks.lat) or not np.array_equal(
        series.lon, masks.lon
    ):
        raise ValueError("series and masks must use the same grid")

    region_masks = masks.masks
    region_masks["total"] = masks.region_id > 0
    region_names = tuple(masks.region_names) + ("total",)
    values_mm: dict[str, np.ndarray] = {}
    contributing_area: dict[str, np.ndarray] = {}
    valid_fraction: dict[str, np.ndarray] = {}

    for name in region_names:
        region_mask = region_masks[name]
        total_area = float(np.sum(masks.cell_area_m2[region_mask]))
        values = np.full(series.months.size, np.nan, dtype=float)
        areas = np.zeros(series.months.size, dtype=float)
        fractions = np.zeros(series.months.size, dtype=float)
        for month_index in range(series.months.size):
            if not series.valid_month[month_index]:
                areas[month_index] = np.nan
                fractions[month_index] = np.nan
                continue
            finite = region_mask & np.isfinite(series.ewh_mm[month_index])
            if not np.any(finite):
                continue
            areas[month_index] = float(np.sum(masks.cell_area_m2[finite]))
            fractions[month_index] = areas[month_index] / total_area if total_area else np.nan
            values[month_index] = float(
                np.sum(series.ewh_mm[month_index][finite] * masks.cell_area_m2[finite])
                / ocean_area_m2
            )
        values_mm[name] = values
        contributing_area[name] = areas
        valid_fraction[name] = fractions

    return RegionalSeries(
        source_id=series.source_id,
        months=series.months,
        region_names=region_names,
        values_mm=values_mm,
        contributing_area_m2=contributing_area,
        valid_cell_fraction=valid_fraction,
        month_status=series.month_status,
        metadata={
            **series.metadata,
            "ocean_area_m2": float(ocean_area_m2),
            "water_density_kg_m3": float(water_density_kg_m3),
            "density_cancels": True,
            "formula": "sum(ewh_mm * cell_area_m2) / ocean_area_m2",
        },
    )
