"""Per-center interannual processing and multi-center ensemble statistics."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Mapping, Sequence

import numpy as np

from .integrate import RegionalSeries


@dataclass(frozen=True)
class ProcessedRegionalSeries:
    """All intermediate temporal states for one center's regional series."""

    source_id: str
    months: np.ndarray
    region_names: tuple[str, ...]
    raw_mm: Mapping[str, np.ndarray]
    deseasoned_mm: Mapping[str, np.ndarray]
    detrended_mm: Mapping[str, np.ndarray]
    smoothed_mm: Mapping[str, np.ndarray]
    month_status: np.ndarray
    metadata: Mapping[str, Any]

    def __post_init__(self) -> None:
        months = np.asarray(self.months, dtype="U7")
        names = tuple(self.region_names)
        for collection_name in (
            "raw_mm",
            "deseasoned_mm",
            "detrended_mm",
            "smoothed_mm",
        ):
            collection = getattr(self, collection_name)
            if set(collection) != set(names):
                raise ValueError(f"{collection_name} must contain every region")
            normalized = {
                name: np.asarray(collection[name], dtype=float) for name in names
            }
            if any(value.shape != (months.size,) for value in normalized.values()):
                raise ValueError(f"{collection_name} arrays must match the month axis")
            object.__setattr__(self, collection_name, normalized)
        status = np.asarray(self.month_status, dtype="U13")
        if status.shape != (months.size,):
            raise ValueError("month_status must match the month axis")
        object.__setattr__(self, "months", months)
        object.__setattr__(self, "region_names", names)
        object.__setattr__(self, "month_status", status)
        object.__setattr__(self, "metadata", dict(self.metadata))


@dataclass(frozen=True)
class RegionalEnsemble:
    """Month-by-month arithmetic ensemble of processed center series."""

    months: np.ndarray
    region_names: tuple[str, ...]
    mean_mm: Mapping[str, np.ndarray]
    sample_std_mm: Mapping[str, np.ndarray]
    minimum_mm: Mapping[str, np.ndarray]
    maximum_mm: Mapping[str, np.ndarray]
    valid_center_count: Mapping[str, np.ndarray]
    center_ids: tuple[str, ...]
    metadata: Mapping[str, Any]

    def __post_init__(self) -> None:
        months = np.asarray(self.months, dtype="U7")
        names = tuple(self.region_names)
        for collection_name, dtype in (
            ("mean_mm", float),
            ("sample_std_mm", float),
            ("minimum_mm", float),
            ("maximum_mm", float),
            ("valid_center_count", int),
        ):
            collection = getattr(self, collection_name)
            if set(collection) != set(names):
                raise ValueError(f"{collection_name} must contain every region")
            normalized = {
                name: np.asarray(collection[name], dtype=dtype) for name in names
            }
            if any(value.shape != (months.size,) for value in normalized.values()):
                raise ValueError(f"{collection_name} arrays must match the month axis")
            object.__setattr__(self, collection_name, normalized)
        object.__setattr__(self, "months", months)
        object.__setattr__(self, "region_names", names)
        object.__setattr__(self, "center_ids", tuple(self.center_ids))
        object.__setattr__(self, "metadata", dict(self.metadata))


def centered_three_month_mean(values: np.ndarray) -> np.ndarray:
    """Strict centered mean; a single non-finite member invalidates the window."""
    values = np.asarray(values, dtype=float)
    if values.ndim != 1:
        raise ValueError("centered smoothing expects a one-dimensional series")
    output = np.full(values.shape, np.nan, dtype=float)
    for index in range(1, values.size - 1):
        window = values[index - 1 : index + 2]
        if np.isfinite(window).all():
            output[index] = float(np.mean(window))
    return output


def _decimal_year(months: np.ndarray) -> np.ndarray:
    years = np.asarray([int(month[:4]) for month in months], dtype=float)
    month_numbers = np.asarray([int(month[5:7]) for month in months], dtype=float)
    return years + (month_numbers - 0.5) / 12.0


def process_interannual(series: RegionalSeries) -> ProcessedRegionalSeries:
    """Deseason, OLS-detrend, then apply a strict centered 3-month mean."""
    month_numbers = np.asarray([int(month[5:7]) for month in series.months])
    decimal_year = _decimal_year(series.months)
    raw: dict[str, np.ndarray] = {}
    deseasoned: dict[str, np.ndarray] = {}
    detrended: dict[str, np.ndarray] = {}
    smoothed: dict[str, np.ndarray] = {}
    slopes: dict[str, float | None] = {}

    for name in series.region_names:
        values = np.asarray(series.values_mm[name], dtype=float).copy()
        seasonal_removed = np.full(values.shape, np.nan, dtype=float)
        for calendar_month in range(1, 13):
            selection = (month_numbers == calendar_month) & np.isfinite(values)
            if np.any(selection):
                seasonal_removed[selection] = values[selection] - float(
                    np.mean(values[selection])
                )

        finite = np.isfinite(seasonal_removed)
        trend_removed = np.full(values.shape, np.nan, dtype=float)
        if np.count_nonzero(finite) >= 2:
            design = np.column_stack(
                (np.ones(np.count_nonzero(finite)), decimal_year[finite])
            )
            coefficients, _, _, _ = np.linalg.lstsq(
                design, seasonal_removed[finite], rcond=None
            )
            trend_removed[finite] = seasonal_removed[finite] - design @ coefficients
            slopes[name] = float(coefficients[1])
        else:
            trend_removed[finite] = seasonal_removed[finite]
            slopes[name] = None

        raw[name] = values
        deseasoned[name] = seasonal_removed
        detrended[name] = trend_removed
        smoothed[name] = centered_three_month_mean(trend_removed)

    metadata = dict(series.metadata)
    metadata.update(
        {
            "processing_order": [
                "monthly_climatology",
                "ols_detrend",
                "centered_3_month_mean",
            ],
            "ols_slope_mm_per_year": slopes,
            "smoothing_requires_three_finite_months": True,
        }
    )
    return ProcessedRegionalSeries(
        source_id=series.source_id,
        months=series.months,
        region_names=series.region_names,
        raw_mm=raw,
        deseasoned_mm=deseasoned,
        detrended_mm=detrended,
        smoothed_mm=smoothed,
        month_status=series.month_status,
        metadata=metadata,
    )


def _expected_event_months(start: str, end: str) -> np.ndarray:
    first = np.datetime64(start, "M")
    last = np.datetime64(end, "M")
    return np.asarray(
        [
            str(item)
            for item in np.arange(
                first,
                last + np.timedelta64(1, "M"),
                dtype="datetime64[M]",
            )
        ]
    )


def combine_centers(
    series_by_center: Mapping[str, ProcessedRegionalSeries],
    *,
    paper_mode: bool = False,
    events: Sequence[Mapping[str, str]] | None = None,
) -> RegionalEnsemble:
    """Average each center's already-processed regional series."""
    if not series_by_center:
        raise ValueError("at least one center is required")
    center_ids = tuple(series_by_center)
    first = series_by_center[center_ids[0]]
    for center_id, series in series_by_center.items():
        if not np.array_equal(series.months, first.months):
            raise ValueError(f"center {center_id} uses a different month axis")
        if series.region_names != first.region_names:
            raise ValueError(f"center {center_id} uses different regional names")

    means: dict[str, np.ndarray] = {}
    standard_deviations: dict[str, np.ndarray] = {}
    minima: dict[str, np.ndarray] = {}
    maxima: dict[str, np.ndarray] = {}
    counts: dict[str, np.ndarray] = {}
    for name in first.region_names:
        stack = np.stack(
            [series_by_center[center].smoothed_mm[name] for center in center_ids]
        )
        finite = np.isfinite(stack)
        count = np.sum(finite, axis=0)
        total = np.sum(np.where(finite, stack, 0.0), axis=0)
        mean = np.full(first.months.size, np.nan, dtype=float)
        np.divide(total, count, out=mean, where=count > 0)
        sample_std = np.full(first.months.size, np.nan, dtype=float)
        for index in np.flatnonzero(count >= 2):
            sample = stack[finite[:, index], index]
            sample_std[index] = float(np.std(sample, ddof=1))
        minimum = np.full(first.months.size, np.nan, dtype=float)
        maximum = np.full(first.months.size, np.nan, dtype=float)
        for index in np.flatnonzero(count > 0):
            sample = stack[finite[:, index], index]
            minimum[index] = float(np.min(sample))
            maximum[index] = float(np.max(sample))
        means[name] = mean
        standard_deviations[name] = sample_std
        minima[name] = minimum
        maxima[name] = maximum
        counts[name] = count

    ensemble = RegionalEnsemble(
        months=first.months,
        region_names=first.region_names,
        mean_mm=means,
        sample_std_mm=standard_deviations,
        minimum_mm=minima,
        maximum_mm=maxima,
        valid_center_count=counts,
        center_ids=center_ids,
        metadata={
            "center_processing_precedes_averaging": True,
            "ensemble_statistic": "arithmetic_mean",
            "spread_statistic": "sample_standard_deviation",
        },
    )

    if paper_mode:
        month_lookup = set(ensemble.months.tolist())
        for event in events or ():
            expected = _expected_event_months(event["start"], event["end"])
            if any(month not in month_lookup for month in expected):
                raise ValueError(
                    f"paper event window contains missing values: {event['id']}"
                )
            selection = np.isin(ensemble.months, expected)
            for name in ensemble.region_names:
                if not np.isfinite(ensemble.mean_mm[name][selection]).all():
                    raise ValueError(
                        f"paper event window contains missing values: "
                        f"{event['id']} {name}"
                    )
    return ensemble

