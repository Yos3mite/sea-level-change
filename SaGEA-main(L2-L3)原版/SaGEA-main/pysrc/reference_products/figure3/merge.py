"""Observed-first merging of monthly Mascon observations and reconstructions."""

from __future__ import annotations

from collections.abc import Iterable
from typing import Any

import numpy as np

from .types import MonthlyGridSeries


def _require_same_grid(first: MonthlyGridSeries, second: MonthlyGridSeries) -> None:
    if not np.array_equal(first.lat, second.lat) or not np.array_equal(
        first.lon, second.lon
    ):
        raise ValueError("monthly grid products must use the same latitude/longitude grid")


def _month_axis(start: str, end: str) -> np.ndarray:
    start_month = np.datetime64(start, "M")
    end_month = np.datetime64(end, "M")
    if end_month < start_month:
        raise ValueError("end month precedes start month")
    return np.asarray(
        [
            str(item)
            for item in np.arange(
                start_month,
                end_month + np.timedelta64(1, "M"),
                dtype="datetime64[M]",
            )
        ]
    )


def align_reconstruction_baseline(
    observed: MonthlyGridSeries,
    reconstructed: MonthlyGridSeries,
    overlap_months: Iterable[str],
    *,
    minimum_overlap_months: int = 12,
) -> tuple[MonthlyGridSeries, dict[str, Any]]:
    """Add the median observed-minus-reconstructed bias at each grid cell."""
    _require_same_grid(observed, reconstructed)
    observed_index = {month: index for index, month in enumerate(observed.months)}
    reconstructed_index = {
        month: index for index, month in enumerate(reconstructed.months)
    }
    requested = list(overlap_months)
    common = [
        month
        for month in requested
        if month in observed_index
        and month in reconstructed_index
        and observed.valid_month[observed_index[month]]
        and reconstructed.valid_month[reconstructed_index[month]]
    ]
    if not common:
        raise ValueError(
            f"baseline alignment needs {minimum_overlap_months} finite overlap months"
        )

    observed_overlap = np.stack(
        [observed.ewh_mm[observed_index[month]] for month in common]
    )
    reconstructed_overlap = np.stack(
        [reconstructed.ewh_mm[reconstructed_index[month]] for month in common]
    )
    differences = observed_overlap - reconstructed_overlap
    finite_count = np.sum(np.isfinite(differences), axis=0)
    supported = finite_count >= minimum_overlap_months
    if not np.any(supported):
        raise ValueError(
            f"no grid cells have {minimum_overlap_months} finite overlap months"
        )

    with np.errstate(all="ignore"):
        bias = np.nanmedian(differences, axis=0)
    bias[~supported] = np.nan
    aligned_values = reconstructed.ewh_mm + bias[np.newaxis, :, :]
    valid_month = reconstructed.valid_month & np.any(
        np.isfinite(aligned_values), axis=(1, 2)
    )
    aligned_values[~valid_month] = np.nan
    month_status = reconstructed.month_status.copy()
    month_status[~valid_month] = "missing"

    before_finite = np.isfinite(differences)
    after_differences = observed_overlap - (
        reconstructed_overlap + bias[np.newaxis, :, :]
    )
    after_finite = np.isfinite(after_differences)
    finite_bias = bias[np.isfinite(bias)]
    diagnostics = {
        "overlap_months": common,
        "minimum_overlap_months": minimum_overlap_months,
        "finite_cell_fraction": float(np.mean(supported)),
        "bias_min_mm": float(np.min(finite_bias)),
        "bias_mean_mm": float(np.mean(finite_bias)),
        "bias_max_mm": float(np.max(finite_bias)),
        "splice_rms_before_mm": float(
            np.sqrt(np.mean(np.square(differences[before_finite])))
        ),
        "splice_rms_after_mm": float(
            np.sqrt(np.mean(np.square(after_differences[after_finite])))
        ),
    }
    metadata = dict(reconstructed.metadata)
    metadata["baseline_alignment"] = diagnostics
    return (
        MonthlyGridSeries(
            source_id=reconstructed.source_id,
            months=reconstructed.months,
            lat=reconstructed.lat,
            lon=reconstructed.lon,
            ewh_mm=aligned_values,
            valid_month=valid_month,
            month_status=month_status,
            metadata=metadata,
        ),
        diagnostics,
    )


def fill_missing_months(
    observed: MonthlyGridSeries,
    reconstructed: MonthlyGridSeries,
    start: str,
    end: str,
    *,
    reconstruction_end: str = "2022-12",
) -> MonthlyGridSeries:
    """Build a complete month axis while never replacing an observed month."""
    _require_same_grid(observed, reconstructed)
    months = _month_axis(start, end)
    observed_index = {month: index for index, month in enumerate(observed.months)}
    reconstructed_index = {
        month: index for index, month in enumerate(reconstructed.months)
    }
    values = np.full(
        (months.size, observed.lat.size, observed.lon.size), np.nan, dtype=float
    )
    valid = np.zeros(months.size, dtype=bool)
    statuses = np.full(months.size, "missing", dtype="U13")
    used_reconstruction: list[str] = []

    for output_index, month in enumerate(months):
        original_index = observed_index.get(month)
        if original_index is not None and observed.valid_month[original_index]:
            values[output_index] = observed.ewh_mm[original_index]
            valid[output_index] = True
            statuses[output_index] = observed.month_status[original_index]
            continue

        replacement_index = reconstructed_index.get(month)
        if replacement_index is None or not reconstructed.valid_month[replacement_index]:
            continue
        if month > reconstruction_end:
            raise ValueError(
                f"cannot use a reconstructed month after registered reconstruction end "
                f"{reconstruction_end}: {month}"
            )
        values[output_index] = reconstructed.ewh_mm[replacement_index]
        valid[output_index] = True
        statuses[output_index] = "reconstructed"
        used_reconstruction.append(month)

    metadata = dict(observed.metadata)
    metadata.update(
        {
            "reconstruction_source_id": reconstructed.source_id,
            "reconstruction_end": reconstruction_end,
            "reconstructed_months": used_reconstruction,
        }
    )
    return MonthlyGridSeries(
        source_id=observed.source_id,
        months=months,
        lat=observed.lat,
        lon=observed.lon,
        ewh_mm=values,
        valid_month=valid,
        month_status=statuses,
        metadata=metadata,
    )
