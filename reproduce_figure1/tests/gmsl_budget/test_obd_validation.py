import numpy as np
import pandas as pd
import pytest

from gmsl_budget.models import MonthlySeries
from gmsl_budget.validation import compare_obd


def _series(months, values, name):
    return MonthlySeries(pd.to_datetime([f"{month}-15" for month in months]), np.array(values), name, "mm", {})


def test_obd_validation_aligns_months_without_interpolation():
    reference = _series(["2020-01", "2020-03"], [1.0, 3.0], "reference")
    computed = _series(["2020-02", "2020-03"], [8.0, 10.0], "computed")
    metrics = compare_obd(reference, computed)
    assert metrics.n_common == 1
    assert metrics.common_months == ("2020-03",)
    assert metrics.mean_offset_mm == pytest.approx(7.0)


def test_obd_validation_removes_only_common_mean_for_shape_metrics():
    reference = _series(["2020-01", "2020-02", "2020-03"], [1.0, 2.0, 4.0], "reference")
    computed = _series(["2020-01", "2020-02", "2020-03"], [8.0, 9.0, 11.0], "computed")
    metrics = compare_obd(reference, computed)
    assert metrics.mean_offset_mm == pytest.approx(7.0)
    assert metrics.rms_difference_mm == pytest.approx(0.0, abs=1e-12)
    assert metrics.correlation == pytest.approx(1.0)
    assert metrics.trend_difference_mm_per_year == pytest.approx(0.0, abs=1e-12)
