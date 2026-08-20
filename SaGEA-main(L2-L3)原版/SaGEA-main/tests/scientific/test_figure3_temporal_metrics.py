import numpy as np
import pytest

from pysrc.reference_products.figure3.integrate import RegionalSeries
from pysrc.reference_products.figure3.metrics import event_metrics
from pysrc.reference_products.figure3.temporal import (
    ProcessedRegionalSeries,
    RegionalEnsemble,
    centered_three_month_mean,
    combine_centers,
    process_interannual,
)


REGIONS = (
    "africa",
    "asia",
    "europe",
    "north_america",
    "south_america",
    "oceania",
    "total",
)


def _months(start: str, count: int) -> np.ndarray:
    first = np.datetime64(start, "M")
    return np.asarray([str(item) for item in first + np.arange(count)])


def _regional(months: np.ndarray, values: np.ndarray, source_id: str = "center") -> RegionalSeries:
    mapping = {name: np.asarray(values, dtype=float).copy() for name in REGIONS}
    ones = {name: np.ones(len(months), dtype=float) for name in REGIONS}
    return RegionalSeries(
        source_id=source_id,
        months=months,
        region_names=REGIONS,
        values_mm=mapping,
        contributing_area_m2=ones,
        valid_cell_fraction=ones,
        month_status=np.full(len(months), "observed"),
        metadata={},
    )


def _processed(
    months: np.ndarray,
    values: np.ndarray,
    source_id: str,
) -> ProcessedRegionalSeries:
    mapping = {name: np.asarray(values, dtype=float).copy() for name in REGIONS}
    return ProcessedRegionalSeries(
        source_id=source_id,
        months=months,
        region_names=REGIONS,
        raw_mm=mapping,
        deseasoned_mm=mapping,
        detrended_mm=mapping,
        smoothed_mm=mapping,
        month_status=np.full(len(months), "observed"),
        metadata={},
    )


def _ensemble(months: np.ndarray, values_by_region: dict[str, np.ndarray]) -> RegionalEnsemble:
    mean = {
        name: np.asarray(values_by_region.get(name, np.zeros(len(months))), dtype=float)
        for name in REGIONS
    }
    zeros = {name: np.zeros(len(months), dtype=float) for name in REGIONS}
    counts = {name: np.full(len(months), 3, dtype=int) for name in REGIONS}
    return RegionalEnsemble(
        months=months,
        region_names=REGIONS,
        mean_mm=mean,
        sample_std_mm=zeros,
        minimum_mm=mean,
        maximum_mm=mean,
        valid_center_count=counts,
        center_ids=("CSR", "JPL", "GSFC"),
        metadata={},
    )


def test_centered_three_month_mean_requires_three_finite_consecutive_values():
    """Catch permissive smoothing that bridges a missing GRACE month."""
    result = centered_three_month_mean(np.asarray([1.0, 2.0, 3.0, np.nan, 5.0]))
    assert np.isnan(result[0])
    assert result[1] == pytest.approx(2.0)
    assert np.isnan(result[2:]).all()


def test_process_interannual_removes_monthly_cycle_and_linear_trend():
    """Catch a changed processing order or a trend left in Figure 3 curves."""
    months = _months("2010-01", 60)
    index = np.arange(60, dtype=float)
    seasonal = 4.0 * np.sin(2.0 * np.pi * (index % 12) / 12.0)
    raw = 0.25 * index + seasonal

    processed = process_interannual(_regional(months, raw))

    month_numbers = np.asarray([int(month[5:7]) for month in months])
    for calendar_month in range(1, 13):
        selection = month_numbers == calendar_month
        assert np.mean(processed.deseasoned_mm["total"][selection]) == pytest.approx(
            0.0, abs=1e-12
        )
    decimal_year = np.asarray(
        [int(month[:4]) + (int(month[5:7]) - 0.5) / 12.0 for month in months]
    )
    fitted_slope = np.polyfit(decimal_year, processed.detrended_mm["total"], 1)[0]
    assert fitted_slope == pytest.approx(0.0, abs=1e-12)
    assert processed.metadata["processing_order"] == [
        "monthly_climatology",
        "ols_detrend",
        "centered_3_month_mean",
    ]


def test_process_interannual_retains_a_known_interannual_pulse():
    """Catch temporal processing that removes the event signal with seasonality."""
    months = _months("2010-01", 60)
    index = np.arange(60, dtype=float)
    raw = 0.1 * index + 2.0 * np.cos(2.0 * np.pi * (index % 12) / 12.0)
    raw[30:33] += 9.0

    processed = process_interannual(_regional(months, raw))

    assert np.nanmax(processed.smoothed_mm["total"]) > 5.0


def test_combine_centers_averages_processed_series_and_uses_sample_spread():
    """Catch averaging raw center grids before each center is processed."""
    months = _months("2020-01", 3)
    centers = {
        "CSR": _processed(months, np.asarray([1.0, 2.0, 3.0]), "CSR"),
        "JPL": _processed(months, np.asarray([3.0, 4.0, 5.0]), "JPL"),
    }

    ensemble = combine_centers(centers)

    assert ensemble.mean_mm["total"].tolist() == [2.0, 3.0, 4.0]
    assert ensemble.sample_std_mm["total"].tolist() == pytest.approx(
        [np.sqrt(2.0)] * 3
    )
    assert ensemble.valid_center_count["total"].tolist() == [2, 2, 2]


def test_combine_centers_paper_mode_rejects_internal_event_gap():
    """Catch publication plots that silently bridge an unresolved event-window gap."""
    months = _months("2020-01", 5)
    values = np.asarray([1.0, 2.0, np.nan, 4.0, 5.0])
    centers = {"CSR": _processed(months, values, "CSR")}

    with pytest.raises(ValueError, match="paper event window contains missing values"):
        combine_centers(
            centers,
            paper_mode=True,
            events=[{"id": "event", "start": "2020-01", "end": "2020-05"}],
        )


def test_event_metrics_preserve_reduction_sign():
    """Catch bar heights computed as unsigned magnitudes before comparison."""
    months = np.asarray(["2014-10", "2015-12"])
    ensemble = _ensemble(
        months,
        {"total": np.asarray([1.0, -5.37])},
    )

    rows = event_metrics(
        ensemble,
        events=[{"id": "event", "start": "2014-10", "end": "2015-12"}],
        paper_references={"event": {"total": -6.37}},
    )

    total = next(row for row in rows if row["region"] == "total")
    assert total["change_mm"] == pytest.approx(-6.37)
    assert total["reduction_magnitude_mm"] == pytest.approx(6.37)
    assert total["paper_reference_mm"] == pytest.approx(-6.37)
    assert total["difference_from_paper_mm"] == pytest.approx(0.0)


def test_event_metrics_leave_unreported_paper_regions_empty():
    """Catch visual digitization being presented as a published numeric reference."""
    months = np.asarray(["2023-05", "2023-12"])
    ensemble = _ensemble(months, {"asia": np.asarray([0.0, -1.0])})

    rows = event_metrics(
        ensemble,
        events=[{"id": "event", "start": "2023-05", "end": "2023-12"}],
        paper_references={"event": {}},
    )

    asia = next(row for row in rows if row["region"] == "asia")
    assert asia["paper_reference_mm"] is None
    assert asia["difference_from_paper_mm"] is None
    assert asia["reference_status"] == "not_reported"
