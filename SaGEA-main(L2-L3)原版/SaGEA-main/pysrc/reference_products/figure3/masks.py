"""Continental masks, geodesic coastal exclusion, and exact grid-cell areas."""

from __future__ import annotations

from dataclasses import dataclass
import json
from pathlib import Path
from typing import Any, Mapping
import warnings

import geopandas as gpd
from netCDF4 import Dataset
import numpy as np
from scipy.spatial import cKDTree
import shapely
from shapely.geometry import MultiPolygon, Polygon


REGION_NAMES = (
    "africa",
    "asia",
    "europe",
    "north_america",
    "south_america",
    "oceania",
)
_NATURAL_EARTH_CONTINENTS = (
    "Africa",
    "Asia",
    "Europe",
    "North America",
    "South America",
    "Oceania",
)


@dataclass(frozen=True)
class ContinentMaskSet:
    """Six mutually exclusive buffered continental masks on one grid."""

    lat: np.ndarray
    lon: np.ndarray
    region_names: tuple[str, ...]
    region_id: np.ndarray
    land_mask: np.ndarray
    distance_to_coast_km: np.ndarray
    coastal_buffer_excluded: np.ndarray
    cell_area_m2: np.ndarray
    metadata: Mapping[str, Any]

    def __post_init__(self) -> None:
        lat = np.asarray(self.lat, dtype=float)
        lon = np.asarray(self.lon, dtype=float)
        shape = (lat.size, lon.size)
        if lat.ndim != 1 or lon.ndim != 1:
            raise ValueError("mask latitude/longitude must be one-dimensional")
        if np.any(np.diff(lat) <= 0) or np.any(np.diff(lon) <= 0):
            raise ValueError("mask latitude/longitude must be strictly increasing")
        arrays = {
            "region_id": np.asarray(self.region_id, dtype=np.int16),
            "land_mask": np.asarray(self.land_mask, dtype=bool),
            "distance_to_coast_km": np.asarray(
                self.distance_to_coast_km, dtype=float
            ),
            "coastal_buffer_excluded": np.asarray(
                self.coastal_buffer_excluded, dtype=bool
            ),
            "cell_area_m2": np.asarray(self.cell_area_m2, dtype=float),
        }
        for name, value in arrays.items():
            if value.shape != shape:
                raise ValueError(f"{name} shape must be {shape}")
        if not self.region_names:
            raise ValueError("region_names must be non-empty")
        if np.any(arrays["region_id"] < 0) or np.any(
            arrays["region_id"] > len(self.region_names)
        ):
            raise ValueError("region_id contains an unknown region")
        if np.any((arrays["region_id"] > 0) & ~arrays["land_mask"]):
            raise ValueError("continental regions must be a subset of land_mask")
        if np.any(arrays["coastal_buffer_excluded"] & ~arrays["land_mask"]):
            raise ValueError("coastal exclusions must be land cells")
        if np.any(arrays["cell_area_m2"] <= 0):
            raise ValueError("cell areas must be positive")

        object.__setattr__(self, "lat", lat)
        object.__setattr__(self, "lon", lon)
        object.__setattr__(self, "region_names", tuple(self.region_names))
        for name, value in arrays.items():
            object.__setattr__(self, name, value)
        object.__setattr__(self, "metadata", dict(self.metadata))

    @property
    def masks(self) -> dict[str, np.ndarray]:
        return {
            name: self.region_id == index
            for index, name in enumerate(self.region_names, start=1)
        }


def _cell_edges(coordinates: np.ndarray, *, clip: tuple[float, float] | None) -> np.ndarray:
    if coordinates.size == 0:
        raise ValueError("grid coordinates must be non-empty")
    if coordinates.size == 1:
        edges = np.asarray([coordinates[0] - 0.5, coordinates[0] + 0.5])
    else:
        midpoints = (coordinates[:-1] + coordinates[1:]) / 2.0
        edges = np.concatenate(
            (
                [coordinates[0] - (midpoints[0] - coordinates[0])],
                midpoints,
                [coordinates[-1] + (coordinates[-1] - midpoints[-1])],
            )
        )
    if clip is not None:
        edges = np.clip(edges, clip[0], clip[1])
    return edges


def cell_areas_m2(
    lat: np.ndarray,
    lon: np.ndarray,
    radius_m: float = 6_371_000.0,
) -> np.ndarray:
    """Return exact spherical areas for cells defined by coordinate midpoints."""
    lat = np.asarray(lat, dtype=float)
    lon = np.asarray(lon, dtype=float)
    if np.any(np.diff(lat) <= 0) or np.any(np.diff(lon) <= 0):
        raise ValueError("latitude and longitude must be strictly increasing")
    lat_edges = np.deg2rad(_cell_edges(lat, clip=(-90.0, 90.0)))
    lon_edges = np.deg2rad(_cell_edges(lon, clip=None))
    latitude_factor = np.sin(lat_edges[1:]) - np.sin(lat_edges[:-1])
    longitude_width = np.diff(lon_edges)
    return radius_m**2 * latitude_factor[:, np.newaxis] * longitude_width


def _resolve_column(frame: gpd.GeoDataFrame, candidates: tuple[str, ...]) -> str:
    lookup = {str(column).lower(): str(column) for column in frame.columns}
    for candidate in candidates:
        if candidate.lower() in lookup:
            return lookup[candidate.lower()]
    raise ValueError(f"none of the required columns exist: {candidates}")


def _polygon_exterior_vertices(geometry: Any, maximum_segment_deg: float) -> np.ndarray:
    densified = shapely.segmentize(geometry, max_segment_length=maximum_segment_deg)
    polygons: list[Polygon]
    if isinstance(densified, Polygon):
        polygons = [densified]
    elif isinstance(densified, MultiPolygon):
        polygons = list(densified.geoms)
    else:
        polygons = [part for part in shapely.get_parts(densified) if isinstance(part, Polygon)]
    coordinates = [np.asarray(polygon.exterior.coords, dtype=float) for polygon in polygons]
    if not coordinates:
        raise ValueError("land geometry has no exterior coastline vertices")
    return np.concatenate(coordinates, axis=0)


def _unit_sphere(longitude_deg: np.ndarray, latitude_deg: np.ndarray) -> np.ndarray:
    longitude = np.deg2rad(longitude_deg)
    latitude = np.deg2rad(latitude_deg)
    cosine_latitude = np.cos(latitude)
    return np.column_stack(
        (
            cosine_latitude * np.cos(longitude),
            cosine_latitude * np.sin(longitude),
            np.sin(latitude),
        )
    )


def _distance_to_exterior_coast_km(
    geometry: Any,
    longitude_deg: np.ndarray,
    latitude_deg: np.ndarray,
    *,
    radius_m: float,
    coastline_densify_deg: float,
) -> np.ndarray:
    vertices = _polygon_exterior_vertices(geometry, coastline_densify_deg)
    tree = cKDTree(_unit_sphere(vertices[:, 0], vertices[:, 1]))
    chord_distance, _ = tree.query(
        _unit_sphere(longitude_deg, latitude_deg), workers=-1
    )
    central_angle = 2.0 * np.arcsin(np.clip(chord_distance / 2.0, 0.0, 1.0))
    return central_angle * radius_m / 1000.0


def build_continent_masks(
    vector_path: Path,
    target_lat: np.ndarray,
    target_lon: np.ndarray,
    *,
    coastal_exclusion_km: float = 300.0,
    coastline_densify_deg: float = 0.25,
    radius_m: float = 6_371_000.0,
) -> ContinentMaskSet:
    """Build the six Jin-style continental domains from Admin-0 polygons."""
    target_lat = np.asarray(target_lat, dtype=float)
    target_lon = np.asarray(target_lon, dtype=float)
    frame = gpd.read_file(vector_path)
    if frame.crs is not None and not frame.crs.is_geographic:
        frame = frame.to_crs("EPSG:4326")
    continent_column = _resolve_column(frame, ("CONTINENT", "continent"))
    admin_column = _resolve_column(
        frame, ("ADMIN", "admin", "NAME", "name", "SOVEREIGNT")
    )
    frame = frame.loc[frame.geometry.notna() & ~frame.geometry.is_empty].copy()
    full_land_geometry = frame.geometry.union_all()

    admin_names = frame[admin_column].fillna("").astype(str).str.casefold()
    continental_source = frame.loc[~admin_names.eq("greenland")].copy()
    geometries: list[Any] = []
    for label in _NATURAL_EARTH_CONTINENTS:
        selected = continental_source.loc[
            continental_source[continent_column].astype(str) == label
        ]
        geometries.append(selected.geometry.union_all() if not selected.empty else None)

    longitude_mesh, latitude_mesh = np.meshgrid(target_lon, target_lat)
    points = shapely.points(longitude_mesh.ravel(), latitude_mesh.ravel())
    land_mask = np.asarray(
        shapely.covers(full_land_geometry, points), dtype=bool
    ).reshape(longitude_mesh.shape)

    raw_masks: list[np.ndarray] = []
    for geometry in geometries:
        if geometry is None or geometry.is_empty:
            raw_masks.append(np.zeros(longitude_mesh.shape, dtype=bool))
        else:
            raw_masks.append(
                np.asarray(shapely.covers(geometry, points), dtype=bool).reshape(
                    longitude_mesh.shape
                )
            )
    membership_count = np.sum(np.stack(raw_masks, axis=0), axis=0)
    if np.any(membership_count > 1):
        raise ValueError("continent masks overlap at one or more grid-cell centers")

    region_id = np.zeros(longitude_mesh.shape, dtype=np.int16)
    for identifier, mask in enumerate(raw_masks, start=1):
        region_id[mask] = identifier

    distance_to_coast_km = np.full(longitude_mesh.shape, np.nan, dtype=float)
    if np.any(land_mask):
        land_flat = land_mask.ravel()
        distances = _distance_to_exterior_coast_km(
            full_land_geometry,
            longitude_mesh.ravel()[land_flat],
            latitude_mesh.ravel()[land_flat],
            radius_m=radius_m,
            coastline_densify_deg=coastline_densify_deg,
        )
        distance_to_coast_km.ravel()[land_flat] = distances

    coastal_buffer_excluded = (
        (region_id > 0) & (distance_to_coast_km < coastal_exclusion_km)
    )
    region_id[coastal_buffer_excluded] = 0
    return ContinentMaskSet(
        lat=target_lat,
        lon=target_lon,
        region_names=REGION_NAMES,
        region_id=region_id,
        land_mask=land_mask,
        distance_to_coast_km=distance_to_coast_km,
        coastal_buffer_excluded=coastal_buffer_excluded,
        cell_area_m2=cell_areas_m2(target_lat, target_lon, radius_m),
        metadata={
            "vector_path": str(vector_path),
            "coastal_exclusion_km": float(coastal_exclusion_km),
            "coastline_densify_deg": float(coastline_densify_deg),
            "radius_m": float(radius_m),
            "greenland_excluded": True,
            "antarctica_excluded": True,
            "coast_definition": "dissolved polygon exterior rings; interior lake rings ignored",
        },
    )


def write_mask_netcdf(
    masks: ContinentMaskSet,
    path: Path,
    *,
    source_hashes: Mapping[str, str] | None = None,
) -> Path:
    """Write the complete auditable mask state to NetCDF."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with Dataset(path, "w", format="NETCDF4") as dataset:
        with warnings.catch_warnings():
            warnings.filterwarnings(
                "ignore",
                message="Setting the shape on a NumPy array has been deprecated",
                category=DeprecationWarning,
            )
            dataset.createDimension("lat", masks.lat.size)
            dataset.createDimension("lon", masks.lon.size)
            dataset.createVariable("lat", "f8", ("lat",))[:] = masks.lat
            dataset.createVariable("lon", "f8", ("lon",))[:] = masks.lon
            dataset.createVariable("region_id", "i2", ("lat", "lon"))[:] = (
                masks.region_id
            )
            dataset.createVariable("land_mask", "i1", ("lat", "lon"))[:] = (
                masks.land_mask
            )
            dataset.createVariable(
                "distance_to_coast_km", "f8", ("lat", "lon"), fill_value=np.nan
            )[:] = masks.distance_to_coast_km
            dataset.createVariable(
                "coastal_buffer_excluded", "i1", ("lat", "lon")
            )[:] = masks.coastal_buffer_excluded
            dataset.createVariable("cell_area_m2", "f8", ("lat", "lon"))[:] = (
                masks.cell_area_m2
            )
        dataset.region_names = json.dumps(
            {str(index): name for index, name in enumerate(masks.region_names, start=1)},
            sort_keys=True,
        )
        dataset.mask_metadata = json.dumps(masks.metadata, sort_keys=True)
        dataset.source_hashes = json.dumps(dict(source_hashes or {}), sort_keys=True)
    return path
