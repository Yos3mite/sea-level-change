import json
from pathlib import Path

from netCDF4 import Dataset
import numpy as np
import pytest

from pysrc.reference_products.figure3.integrate import (
    RegionalSeries,
    integrate_regions,
)
from pysrc.reference_products.figure3.masks import (
    ContinentMaskSet,
    build_continent_masks,
    cell_areas_m2,
    write_mask_netcdf,
)
from pysrc.reference_products.figure3.types import MonthlyGridSeries


def _grid_series(field_mm: np.ndarray) -> MonthlyGridSeries:
    field_mm = np.asarray(field_mm, dtype=float)
    return MonthlyGridSeries(
        source_id="center",
        months=np.asarray(["2020-01"] * field_mm.shape[0]),
        lat=np.arange(field_mm.shape[1], dtype=float),
        lon=np.arange(field_mm.shape[2], dtype=float),
        ewh_mm=field_mm,
        valid_month=np.ones(field_mm.shape[0], dtype=bool),
        month_status=np.full(field_mm.shape[0], "observed"),
        metadata={},
    )


def _mask_set(region_id: np.ndarray, cell_area: np.ndarray) -> ContinentMaskSet:
    region_id = np.asarray(region_id, dtype=np.int16)
    cell_area = np.asarray(cell_area, dtype=float)
    return ContinentMaskSet(
        lat=np.arange(region_id.shape[0], dtype=float),
        lon=np.arange(region_id.shape[1], dtype=float),
        region_names=(
            "africa",
            "asia",
            "europe",
            "north_america",
            "south_america",
            "oceania",
        ),
        region_id=region_id,
        land_mask=region_id > 0,
        distance_to_coast_km=np.full(region_id.shape, 500.0),
        coastal_buffer_excluded=np.zeros(region_id.shape, dtype=bool),
        cell_area_m2=cell_area,
        metadata={},
    )


def _write_geojson(path: Path, features: list[dict]) -> Path:
    path.write_text(
        json.dumps({"type": "FeatureCollection", "features": features}),
        encoding="utf-8",
    )
    return path


def _feature(continent: str, admin: str, coordinates: list[list[float]]) -> dict:
    return {
        "type": "Feature",
        "properties": {"CONTINENT": continent, "ADMIN": admin},
        "geometry": {"type": "Polygon", "coordinates": [coordinates]},
    }


def test_cell_areas_uses_exact_spherical_latitude_edges():
    """Catch cosine-center approximations replacing exact spherical bands."""
    areas = cell_areas_m2(
        lat=np.asarray([-0.5, 0.5]),
        lon=np.asarray([-0.5, 0.5]),
        radius_m=1.0,
    )
    expected = np.deg2rad(1.0) * (np.sin(np.deg2rad(0.0)) - np.sin(np.deg2rad(-1.0)))
    assert areas.shape == (2, 2)
    assert areas[0, 0] == pytest.approx(expected)
    assert areas[0, 0] == pytest.approx(areas[1, 1])


def test_integrate_regions_uses_mass_not_regional_mean():
    """Catch regional-mean EWH being mislabeled as global ESL contribution."""
    masks = _mask_set(region_id=np.asarray([[1, 1]]), cell_area=np.asarray([[1.0, 2.0]]))
    series = _grid_series(np.asarray([[[10.0, 10.0]]]))

    result = integrate_regions(series, masks, ocean_area_m2=6.0)

    assert isinstance(result, RegionalSeries)
    assert result.values_mm["africa"][0] == pytest.approx(5.0)
    assert result.values_mm["total"][0] == pytest.approx(5.0)


def test_total_integration_equals_union_of_six_regions():
    """Catch the Total series using a different land domain from its components."""
    masks = _mask_set(
        region_id=np.asarray([[1, 2, 3, 4, 5, 6]]),
        cell_area=np.ones((1, 6)),
    )
    series = _grid_series(np.asarray([[[1.0, 2.0, 3.0, 4.0, 5.0, 6.0]]]))

    result = integrate_regions(series, masks, ocean_area_m2=2.0)

    components = sum(result.values_mm[name][0] for name in masks.region_names)
    assert result.values_mm["total"][0] == pytest.approx(components)


def test_build_continent_masks_excludes_greenland_and_antarctica(tmp_path: Path):
    """Catch polar ice-sheet mass entering the six-continent TWS total."""
    square = lambda x: [[x - 1, -1], [x + 1, -1], [x + 1, 1], [x - 1, 1], [x - 1, -1]]
    source = _write_geojson(
        tmp_path / "continents.geojson",
        [
            _feature("Africa", "Example Africa", square(0)),
            _feature("North America", "Greenland", square(10)),
            _feature("Antarctica", "Antarctica", square(20)),
        ],
    )

    masks = build_continent_masks(
        source,
        target_lat=np.asarray([0.0]),
        target_lon=np.asarray([0.0, 10.0, 20.0]),
        coastal_exclusion_km=0.0,
    )

    assert masks.land_mask.tolist() == [[True, True, True]]
    assert masks.region_id.tolist() == [[1, 0, 0]]


def test_build_continent_masks_applies_land_side_300km_coastal_exclusion(
    tmp_path: Path,
):
    """Catch a sea-side buffer or a degree-based substitute for 300 km distance."""
    source = _write_geojson(
        tmp_path / "coast.geojson",
        [
            _feature(
                "Africa",
                "Synthetic Africa",
                [[-10, -10], [10, -10], [10, 10], [-10, 10], [-10, -10]],
            )
        ],
    )

    masks = build_continent_masks(
        source,
        target_lat=np.asarray([0.0]),
        target_lon=np.asarray([7.0, 9.0]),
        coastal_exclusion_km=300.0,
        coastline_densify_deg=0.1,
    )

    assert masks.distance_to_coast_km[0, 0] > 300.0
    assert masks.distance_to_coast_km[0, 1] < 300.0
    assert masks.coastal_buffer_excluded.tolist() == [[False, True]]
    assert masks.region_id.tolist() == [[1, 0]]


def test_build_continent_masks_rejects_overlapping_continent_cells(tmp_path: Path):
    """Catch ambiguous double-counting when source continent polygons overlap."""
    square = [[-2, -2], [2, -2], [2, 2], [-2, 2], [-2, -2]]
    source = _write_geojson(
        tmp_path / "overlap.geojson",
        [
            _feature("Africa", "A", square),
            _feature("Asia", "B", square),
        ],
    )

    with pytest.raises(ValueError, match="continent masks overlap"):
        build_continent_masks(
            source,
            target_lat=np.asarray([0.0]),
            target_lon=np.asarray([0.0]),
            coastal_exclusion_km=0.0,
        )


def test_write_mask_netcdf_preserves_auditable_mask_contract(tmp_path: Path):
    """Catch mask outputs that cannot reproduce region membership and exclusions."""
    masks = _mask_set(np.asarray([[1, 0]]), np.asarray([[1.0, 2.0]]))
    path = tmp_path / "masks.nc"

    write_mask_netcdf(masks, path, source_hashes={"natural_earth": "abc"})

    with Dataset(path) as dataset:
        assert set(
            [
                "region_id",
                "land_mask",
                "distance_to_coast_km",
                "coastal_buffer_excluded",
                "cell_area_m2",
            ]
        ).issubset(dataset.variables)
        assert json.loads(dataset.region_names)["1"] == "africa"
        assert json.loads(dataset.source_hashes) == {"natural_earth": "abc"}
