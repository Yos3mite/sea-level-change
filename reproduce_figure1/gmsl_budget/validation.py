from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .models import MonthlySeries
from .trend import decimal_year


@dataclass(frozen=True)
class ValidationMetrics:
    n_common: int
    common_months: tuple[str, ...]
    mean_offset_mm: float
    trend_difference_mm_per_year: float
    correlation: float
    rms_difference_mm: float


def compare_obd(reference: MonthlySeries, computed: MonthlySeries) -> ValidationMetrics:
    if reference.units != "mm" or computed.units != "mm":
        raise ValueError("OBD validation series must use millimetres")
    reference_index = {str(month): index for index, month in enumerate(reference.time.to_period("M"))}
    computed_index = {str(month): index for index, month in enumerate(computed.time.to_period("M"))}
    common = sorted(set(reference_index) & set(computed_index))
    if not common:
        raise ValueError("OBD validation has no common months")
    reference_values = np.array([reference.values[reference_index[month]] for month in common], dtype=np.float64)
    computed_values = np.array([computed.values[computed_index[month]] for month in common], dtype=np.float64)
    valid = np.isfinite(reference_values) & np.isfinite(computed_values)
    common = [month for month, keep in zip(common, valid) if keep]
    reference_values = reference_values[valid]
    computed_values = computed_values[valid]
    if not common:
        raise ValueError("OBD validation has no finite common months")
    mean_offset = float(np.mean(computed_values - reference_values))
    reference_anomaly = reference_values - np.mean(reference_values)
    computed_anomaly = computed_values - np.mean(computed_values)
    rms = float(np.sqrt(np.mean((computed_anomaly - reference_anomaly) ** 2)))
    if len(common) >= 2 and np.std(reference_anomaly) > 0 and np.std(computed_anomaly) > 0:
        correlation = float(np.corrcoef(reference_anomaly, computed_anomaly)[0, 1])
    else:
        correlation = float("nan")
    if len(common) >= 2:
        dates = np.array([np.datetime64(f"{month}-15") for month in common])
        year = decimal_year(dates)
        reference_trend = float(np.polyfit(year - year.mean(), reference_values, 1)[0])
        computed_trend = float(np.polyfit(year - year.mean(), computed_values, 1)[0])
        trend_difference = computed_trend - reference_trend
    else:
        trend_difference = float("nan")
    return ValidationMetrics(
        n_common=len(common),
        common_months=tuple(common),
        mean_offset_mm=mean_offset,
        trend_difference_mm_per_year=trend_difference,
        correlation=correlation,
        rms_difference_mm=rms,
    )
