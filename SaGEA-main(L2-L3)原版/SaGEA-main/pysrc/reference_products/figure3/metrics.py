"""Signed event endpoint metrics for Figure 3 curves and bars."""

from __future__ import annotations

from typing import Any, Mapping, Sequence

import numpy as np

from .temporal import RegionalEnsemble


def event_metrics(
    ensemble: RegionalEnsemble,
    events: Sequence[Mapping[str, str]],
    paper_references: Mapping[str, Mapping[str, float]],
) -> list[dict[str, Any]]:
    """Compute signed endpoint changes without tuning to paper reference values."""
    month_index = {month: index for index, month in enumerate(ensemble.months)}
    rows: list[dict[str, Any]] = []
    for event in events:
        event_id = event["id"]
        start = event["start"]
        end = event["end"]
        if start not in month_index or end not in month_index:
            raise ValueError(f"event endpoints are absent from month axis: {event_id}")
        start_index = month_index[start]
        end_index = month_index[end]
        event_references = paper_references.get(event_id, {})
        for name in ensemble.region_names:
            start_value = float(ensemble.mean_mm[name][start_index])
            end_value = float(ensemble.mean_mm[name][end_index])
            change = end_value - start_value
            reference = event_references.get(name)
            difference = None if reference is None else change - float(reference)
            rows.append(
                {
                    "event_id": event_id,
                    "region": name,
                    "start_month": start,
                    "end_month": end,
                    "start_mm": start_value,
                    "end_mm": end_value,
                    "change_mm": change,
                    "reduction_magnitude_mm": max(-change, 0.0),
                    "paper_reference_mm": None
                    if reference is None
                    else float(reference),
                    "difference_from_paper_mm": difference,
                    "reference_status": "not_reported"
                    if reference is None
                    else "reported",
                    "basis": "ensemble_mean_smoothed_endpoint_difference",
                    "valid_centers_start": int(
                        ensemble.valid_center_count[name][start_index]
                    ),
                    "valid_centers_end": int(
                        ensemble.valid_center_count[name][end_index]
                    ),
                }
            )
    return rows
