from __future__ import annotations

import numpy as np
import pandas as pd

from .grace import ensemble_mean
from .models import MonthlySeries
from .trend import decimal_year, fit_trend


def _diagonal_average(matrix: np.ndarray) -> np.ndarray:
    rows, columns = matrix.shape
    values = np.zeros(rows + columns - 1, dtype=np.float64)
    counts = np.zeros(rows + columns - 1, dtype=np.float64)
    for row in range(rows):
        indices = row + np.arange(columns)
        values[indices] += matrix[row]
        counts[indices] += 1.0
    return values / counts


def iterative_ssa_fill(
    series: MonthlySeries,
    start_month: str,
    end_month: str,
    window: int = 36,
    rank: int = 8,
    max_iterations: int = 200,
    tolerance: float = 1.0e-8,
) -> MonthlySeries:
    if series.units != "mm":
        raise ValueError("SSA input must be expressed in millimetres (mm)")
    periods = pd.period_range(start_month, end_month, freq="M")
    if len(periods) < 4:
        raise ValueError("SSA requires at least four target months")
    if not 2 <= window < len(periods):
        raise ValueError("SSA window must be in [2, number of months)")
    columns = len(periods) - window + 1
    if not 1 <= rank <= min(window, columns):
        raise ValueError("SSA rank exceeds the trajectory matrix dimensions")
    if max_iterations <= 0 or tolerance <= 0:
        raise ValueError("SSA iteration controls must be positive")
    source = {month: value for month, value in zip(series.time.to_period("M"), series.values)}
    observed_values = np.asarray([source.get(month, np.nan) for month in periods], dtype=np.float64)
    observed = np.isfinite(observed_values)
    if int(observed.sum()) < max(window, rank + 1):
        raise ValueError("SSA has too few observed months for the requested window and rank")
    estimate = pd.Series(observed_values).interpolate(limit_direction="both").to_numpy(dtype=np.float64)
    missing = ~observed
    converged = not np.any(missing)
    iterations = 0
    for iterations in range(1, max_iterations + 1):
        trajectory = np.lib.stride_tricks.sliding_window_view(estimate, window).T
        left, singular, right = np.linalg.svd(trajectory, full_matrices=False)
        reconstructed_matrix = (left[:, :rank] * singular[:rank]) @ right[:rank]
        reconstructed = _diagonal_average(reconstructed_matrix)
        updated = estimate.copy()
        updated[missing] = reconstructed[missing]
        updated[observed] = observed_values[observed]
        denominator = max(float(np.linalg.norm(estimate[missing])), 1.0)
        relative_change = float(np.linalg.norm(updated[missing] - estimate[missing]) / denominator)
        estimate = updated
        if relative_change <= tolerance:
            converged = True
            break
    if not converged:
        raise RuntimeError(f"SSA did not converge in {max_iterations} iterations")
    metadata = dict(series.metadata)
    metadata.update(
        {
            "gap_fill_method": "iterative singular spectrum analysis",
            "gap_fill_scope": "missing months only; observed values preserved",
            "gap_fill_source_status": "time-series fallback, not the published Xie and Yi gridded product",
            "ssa_window_months": int(window),
            "ssa_rank": int(rank),
            "ssa_iterations": int(iterations),
            "ssa_tolerance": float(tolerance),
            "filled_month_count": int(missing.sum()),
            "filled_months": periods[missing].astype(str).tolist(),
        }
    )
    return MonthlySeries(
        periods.to_timestamp() + pd.Timedelta(days=14),
        estimate,
        f"{series.name}_ssa_filled",
        "mm",
        metadata,
    )


def build_ssa_ensemble(
    series_by_center: dict[str, MonthlySeries],
    start_month: str,
    end_month: str,
    window: int = 36,
    rank: int = 8,
    max_iterations: int = 1000,
    tolerance: float = 1.0e-5,
) -> tuple[dict[str, MonthlySeries], MonthlySeries]:
    """Gap-fill each processing centre independently, then average centers."""
    if set(series_by_center) != {"CSR", "JPL", "GSFC"}:
        raise ValueError("the Jin ensemble requires exactly CSR, JPL, and GSFC")
    filled = {
        center: iterative_ssa_fill(
            series_by_center[center],
            start_month,
            end_month,
            window=window,
            rank=rank,
            max_iterations=max_iterations,
            tolerance=tolerance,
        )
        for center in sorted(series_by_center)
    }
    base = ensemble_mean(filled)
    metadata = dict(base.metadata)
    metadata.update(
        {
            "processing_order": "SSA separately by center, then arithmetic mean",
            "jin_domain": "global ocean; no coastal buffer",
            "gia_policy": "product GIA correction retained; no second GIA correction",
        }
    )
    return filled, MonthlySeries(base.time, base.values, base.name, base.units, metadata)


def prepare_figure1_series(
    series: MonthlySeries,
    remove_trend: bool,
    running_window: int | None,
    hac_lags: int = 12,
) -> MonthlySeries:
    result = fit_trend(series, hac_lags=hac_lags)
    year = decimal_year(series.time)
    seasonal = (
        result.annual_sin_mm * np.sin(2.0 * np.pi * year)
        + result.annual_cos_mm * np.cos(2.0 * np.pi * year)
        + result.semiannual_sin_mm * np.sin(4.0 * np.pi * year)
        + result.semiannual_cos_mm * np.cos(4.0 * np.pi * year)
    )
    values = series.values - seasonal
    if remove_trend:
        valid = np.isfinite(series.values)
        centered_year = year - float(np.mean(year[valid]))
        values = values - result.trend_mm_per_year * centered_year
    values = values - float(np.nanmean(values))
    if running_window is not None:
        if running_window < 1 or running_window % 2 == 0:
            raise ValueError("running_window must be a positive odd number")
        values = (
            pd.Series(values)
            .rolling(window=running_window, center=True, min_periods=running_window)
            .mean()
            .to_numpy(dtype=np.float64)
        )
    metadata = dict(series.metadata)
    metadata.update(
        {
            "seasonal_components_removed": True,
            "linear_trend_removed": bool(remove_trend),
            "running_mean_months": running_window,
            "figure_vertical_reference": "finite-series mean removed",
            "parent_series": series.name,
        }
    )
    suffix = "deseasoned_detrended" if remove_trend else "deseasoned"
    return MonthlySeries(series.time, values, f"{series.name}_{suffix}", "mm", metadata)
