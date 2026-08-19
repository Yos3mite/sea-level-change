import numpy as np
import pandas as pd
import pytest

from gmsl_budget.figure1_mass import build_ssa_ensemble, iterative_ssa_fill, prepare_figure1_series
from gmsl_budget.models import MonthlySeries
from gmsl_budget.trend import decimal_year, fit_trend


def _series(time, values, name="ocean_mass_csr"):
    return MonthlySeries(pd.DatetimeIndex(time), np.asarray(values, dtype=np.float64), name, "mm", {})


def test_iterative_ssa_fills_only_missing_months_and_preserves_observations():
    full_time = pd.date_range("2014-01-01", periods=72, freq="MS") + pd.Timedelta(days=14)
    t = np.arange(72, dtype=np.float64)
    truth = 0.1 * t + 2.0 * np.sin(2.0 * np.pi * t / 12.0)
    observed = np.ones(72, dtype=bool)
    observed[30:41] = False
    source = _series(full_time[observed], truth[observed])
    filled = iterative_ssa_fill(source, "2014-01", "2019-12", window=24, rank=4)
    source_lookup = dict(zip(source.time.to_period("M"), source.values))
    filled_lookup = dict(zip(filled.time.to_period("M"), filled.values))
    assert len(filled.time) == 72
    assert np.all(np.isfinite(filled.values))
    assert [filled_lookup[month] for month in source_lookup] == pytest.approx(list(source_lookup.values()))
    assert np.sqrt(np.mean((filled.values[30:41] - truth[30:41]) ** 2)) < 0.75
    assert filled.metadata["filled_month_count"] == 11
    assert filled.metadata["gap_fill_method"] == "iterative singular spectrum analysis"


def test_prepare_figure1_panel_a_preserves_trend_and_panel_b_removes_it():
    time = pd.date_range("2013-11-01", periods=132, freq="MS") + pd.Timedelta(days=14)
    year = decimal_year(time)
    centered = year - year.mean()
    values = 1.4 * centered + 3.0 * np.sin(2.0 * np.pi * year) + 1.0 * np.cos(4.0 * np.pi * year)
    source = _series(time, values, "ocean_mass_ensemble")
    panel_a = prepare_figure1_series(source, remove_trend=False, running_window=None)
    panel_b = prepare_figure1_series(source, remove_trend=True, running_window=3)
    assert fit_trend(panel_a).trend_mm_per_year == pytest.approx(1.4, abs=1e-10)
    assert fit_trend(panel_b).trend_mm_per_year == pytest.approx(0.0, abs=1e-10)
    assert np.isnan(panel_b.values[0]) and np.isnan(panel_b.values[-1])
    assert panel_a.metadata["seasonal_components_removed"] is True
    assert panel_b.metadata["linear_trend_removed"] is True
    assert panel_b.metadata["running_mean_months"] == 3


def test_build_ssa_ensemble_fills_each_center_before_three_center_mean():
    time = pd.period_range("2018-01", "2020-12", freq="M").to_timestamp() + pd.Timedelta(days=14)
    base = np.sin(np.arange(len(time)) * 2.0 * np.pi / 12.0)
    members = {}
    for center, offset, missing in (("CSR", 0.0, 5), ("JPL", 1.0, 10), ("GSFC", 2.0, 15)):
        values = base + offset
        values[missing] = np.nan
        members[center] = MonthlySeries(
            time,
            values,
            f"ocean_mass_{center.lower()}",
            "mm",
            {"center": center, "mask_hash": "global-mask"},
        )
    filled, ensemble = build_ssa_ensemble(
        members, "2018-01", "2020-12", window=12, rank=4, max_iterations=500, tolerance=1.0e-6
    )
    assert set(filled) == {"CSR", "JPL", "GSFC"}
    assert all(np.isfinite(item.values).all() for item in filled.values())
    assert ensemble.metadata["processing_order"] == "SSA separately by center, then arithmetic mean"
    assert ensemble.values.tolist() == pytest.approx(
        np.mean(np.vstack([filled[center].values for center in sorted(filled)]), axis=0)
    )
