from __future__ import annotations

import numpy as np
import pandas as pd
import statsmodels.api as sm

from .models import MonthlySeries, TrendResult


def decimal_year(time: pd.DatetimeIndex) -> np.ndarray:
    timestamps = pd.DatetimeIndex(pd.to_datetime(time))
    starts = pd.to_datetime(timestamps.year.astype(str) + "-01-01")
    ends = pd.to_datetime((timestamps.year + 1).astype(str) + "-01-01")
    elapsed = (timestamps - starts).total_seconds().to_numpy(dtype=np.float64)
    duration = (ends - starts).total_seconds().to_numpy(dtype=np.float64)
    return timestamps.year.to_numpy(dtype=np.float64) + elapsed / duration


def _missing_months(series: MonthlySeries) -> tuple[str, ...]:
    if len(series.time) == 0:
        return ()
    present = series.time.to_period("M")
    complete = pd.period_range(present.min(), present.max(), freq="M")
    absent = set(complete.astype(str)) - set(present.astype(str))
    invalid = set(present[~np.isfinite(series.values)].astype(str))
    return tuple(sorted(absent | invalid))


def fit_trend(series: MonthlySeries, hac_lags: int = 12) -> TrendResult:
    if series.units != "mm":
        raise ValueError("trend input must be expressed in millimetres (mm)")
    if hac_lags < 0:
        raise ValueError("hac_lags must be non-negative")
    valid = np.isfinite(series.values)
    if int(valid.sum()) < 36:
        raise ValueError("trend requires at least 36 valid months spanning 3 years")
    time = series.time[valid]
    t = decimal_year(time)
    if float(np.ptp(t)) < 3.0:
        raise ValueError("trend requires at least 36 valid months spanning 3 years")
    centered = t - t.mean()
    design = np.column_stack(
        (
            np.ones(len(t), dtype=np.float64),
            centered,
            np.sin(2.0 * np.pi * t),
            np.cos(2.0 * np.pi * t),
            np.sin(4.0 * np.pi * t),
            np.cos(4.0 * np.pi * t),
        )
    )
    if np.linalg.matrix_rank(design) != design.shape[1]:
        raise ValueError("trend design matrix is rank deficient")
    condition_number = float(np.linalg.cond(design))
    if not np.isfinite(condition_number) or condition_number > 1.0e12:
        raise ValueError(f"trend design matrix is ill-conditioned: {condition_number:.3e}")
    model = sm.OLS(series.values[valid], design).fit()
    robust = model.get_robustcov_results(cov_type="HAC", maxlags=hac_lags)
    residuals = np.asarray(model.resid, dtype=np.float64)
    if len(residuals) > 1 and np.std(residuals[:-1]) > 0 and np.std(residuals[1:]) > 0:
        lag1 = float(np.corrcoef(residuals[:-1], residuals[1:])[0, 1])
    else:
        lag1 = float("nan")
    params = np.asarray(model.params, dtype=np.float64)
    return TrendResult(
        series_name=series.name,
        trend_mm_per_year=float(params[1]),
        ols_standard_error=float(model.bse[1]),
        hac_standard_error=float(robust.bse[1]),
        intercept_mm=float(params[0]),
        annual_sin_mm=float(params[2]),
        annual_cos_mm=float(params[3]),
        semiannual_sin_mm=float(params[4]),
        semiannual_cos_mm=float(params[5]),
        residual_lag1_correlation=lag1,
        n_obs=int(valid.sum()),
        start_time=pd.Timestamp(time.min()),
        end_time=pd.Timestamp(time.max()),
        missing_months=_missing_months(series),
    )
