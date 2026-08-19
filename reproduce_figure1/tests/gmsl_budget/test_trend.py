import numpy as np
import pandas as pd
import pytest

from gmsl_budget.models import MonthlySeries
from gmsl_budget.trend import decimal_year, fit_trend


def _independent_decimal_year(time: pd.DatetimeIndex) -> np.ndarray:
    result = []
    for timestamp in time:
        start = pd.Timestamp(timestamp.year, 1, 1)
        end = pd.Timestamp(timestamp.year + 1, 1, 1)
        result.append(timestamp.year + (timestamp - start) / (end - start))
    return np.asarray(result, dtype=np.float64)


def _month_midpoints(start: str, periods: int) -> pd.DatetimeIndex:
    return pd.period_range(start, periods=periods, freq="M").to_timestamp() + pd.Timedelta(days=14)


def test_decimal_year_respects_leap_year_length():
    time = pd.to_datetime(["2019-07-02", "2020-07-02"])
    got = decimal_year(time)
    expected = _independent_decimal_year(time)
    assert got == pytest.approx(expected, abs=1e-14)
    assert got[0] != got[1] - 1.0


def test_fit_trend_recovers_calendar_trend_and_harmonics():
    time = _month_midpoints("2010-01", 120)
    year = _independent_decimal_year(time)
    centered = year - year.mean()
    values = (
        7.0
        + 3.2 * centered
        + 4.0 * np.sin(2 * np.pi * year)
        + 2.5 * np.cos(2 * np.pi * year)
        - 1.5 * np.sin(4 * np.pi * year)
        + 0.75 * np.cos(4 * np.pi * year)
    )
    result = fit_trend(MonthlySeries(time, values, "synthetic", "mm", {}), hac_lags=12)
    assert result.trend_mm_per_year == pytest.approx(3.2, abs=1e-10)
    assert result.annual_sin_mm == pytest.approx(4.0, abs=1e-10)
    assert result.annual_cos_mm == pytest.approx(2.5, abs=1e-10)
    assert result.semiannual_sin_mm == pytest.approx(-1.5, abs=1e-10)
    assert result.semiannual_cos_mm == pytest.approx(0.75, abs=1e-10)
    assert result.n_obs == 120
    assert result.missing_months == ()


@pytest.mark.parametrize("periods", [24, 35])
def test_fit_trend_rejects_short_series(periods):
    series = MonthlySeries(
        _month_midpoints("2020-01", periods),
        np.ones(periods),
        "x",
        "mm",
        {},
    )
    with pytest.raises(ValueError, match="36 valid months"):
        fit_trend(series)


def test_fit_trend_records_missing_months_without_interpolation():
    time = _month_midpoints("2010-01", 60)
    values = np.arange(60, dtype=np.float64)
    values[[4, 19]] = np.nan
    result = fit_trend(MonthlySeries(time, values, "gappy", "mm", {}))
    assert result.n_obs == 58
    assert result.missing_months == ("2010-05", "2011-08")


def test_fit_trend_rejects_non_mm_units():
    series = MonthlySeries(_month_midpoints("2010-01", 60), np.ones(60), "x", "m", {})
    with pytest.raises(ValueError, match="millimetres"):
        fit_trend(series)
