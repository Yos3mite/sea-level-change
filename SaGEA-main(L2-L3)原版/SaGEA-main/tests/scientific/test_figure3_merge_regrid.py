import numpy as np
import pytest

from pysrc.reference_products.figure3.merge import (
    align_reconstruction_baseline,
    fill_missing_months,
)
from pysrc.reference_products.figure3.regrid import nearest_regrid
from pysrc.reference_products.figure3.types import MonthlyGridSeries


def _series(
    months: list[str],
    values: list[float],
    *,
    status: str = "observed",
    source_id: str = "center",
    lat: np.ndarray | None = None,
    lon: np.ndarray | None = None,
) -> MonthlyGridSeries:
    lat = np.asarray([-0.5, 0.5] if lat is None else lat, dtype=float)
    lon = np.asarray([-0.5, 0.5] if lon is None else lon, dtype=float)
    fields = np.empty((len(months), len(lat), len(lon)), dtype=float)
    valid = np.empty(len(months), dtype=bool)
    statuses = np.empty(len(months), dtype="U13")
    for index, value in enumerate(values):
        if np.isnan(value):
            fields[index] = np.nan
            valid[index] = False
            statuses[index] = "missing"
        else:
            fields[index] = value
            valid[index] = True
            statuses[index] = status
    return MonthlyGridSeries(
        source_id=source_id,
        months=np.asarray(months),
        lat=lat,
        lon=lon,
        ewh_mm=fields,
        valid_month=valid,
        month_status=statuses,
        metadata={},
    )


def test_fill_missing_months_never_overwrites_observed_values():
    """Catch a reconstruction replacing a real observed monthly field."""
    observed = _series(["2017-06", "2017-07"], [1.0, np.nan])
    reconstructed = _series(
        ["2017-06", "2017-07"],
        [101.0, 102.0],
        status="reconstructed",
    )

    merged = fill_missing_months(observed, reconstructed, "2017-06", "2017-07")

    assert merged.ewh_mm[:, 0, 0].tolist() == [1.0, 102.0]
    assert merged.month_status.tolist() == ["observed", "reconstructed"]
    assert merged.metadata["reconstructed_months"] == ["2017-07"]


def test_fill_missing_months_preserves_unresolved_gaps():
    """Catch absent reconstruction months being hidden or interpolated."""
    observed = _series(["2017-06"], [1.0])
    reconstructed = _series(["2017-08"], [3.0], status="reconstructed")

    merged = fill_missing_months(observed, reconstructed, "2017-06", "2017-08")

    assert merged.months.tolist() == ["2017-06", "2017-07", "2017-08"]
    assert merged.month_status.tolist() == ["observed", "missing", "reconstructed"]
    assert np.isnan(merged.ewh_mm[1]).all()


def test_fill_missing_months_rejects_reconstruction_after_registered_end():
    """Catch the Xie-Yi product being extrapolated beyond its 2022-12 coverage."""
    observed = _series(["2023-01"], [np.nan])
    reconstructed = _series(["2023-01"], [9.0], status="reconstructed")

    with pytest.raises(ValueError, match="after registered reconstruction end"):
        fill_missing_months(observed, reconstructed, "2023-01", "2023-01")


def test_align_reconstruction_baseline_removes_additive_grid_bias():
    """Catch splices caused by differing product anomaly reference periods."""
    months = [f"2020-{month:02d}" for month in range(1, 13)]
    observed = _series(months, [float(i + 5) for i in range(12)])
    reconstructed = _series(
        months,
        [float(i) for i in range(12)],
        status="reconstructed",
    )

    aligned, diagnostics = align_reconstruction_baseline(
        observed,
        reconstructed,
        overlap_months=months,
    )

    assert np.allclose(aligned.ewh_mm, reconstructed.ewh_mm + 5.0)
    assert diagnostics["bias_min_mm"] == pytest.approx(5.0)
    assert diagnostics["bias_mean_mm"] == pytest.approx(5.0)
    assert diagnostics["bias_max_mm"] == pytest.approx(5.0)
    assert diagnostics["splice_rms_after_mm"] == pytest.approx(0.0)


def test_align_reconstruction_baseline_rejects_no_supported_grid_cells():
    """Catch a baseline correction inferred from fewer than 12 overlap months."""
    months = [f"2020-{month:02d}" for month in range(1, 12)]
    observed = _series(months, [5.0] * 11)
    reconstructed = _series(months, [0.0] * 11, status="reconstructed")

    with pytest.raises(ValueError, match="12 finite overlap months"):
        align_reconstruction_baseline(observed, reconstructed, months)


def test_nearest_regrid_is_cyclic_at_the_dateline():
    """Catch a false blank or wrong neighbor at the longitude seam."""
    source = _series(
        ["2020-01"],
        [0.0],
        lat=np.asarray([-0.5, 0.5]),
        lon=np.asarray([-179.5, 179.5]),
    )
    source.ewh_mm[0, :, 0] = 10.0
    source.ewh_mm[0, :, 1] = 20.0

    target = nearest_regrid(
        source,
        target_lat=np.asarray([-0.5, 0.5]),
        target_lon=np.asarray([-179.6, 179.6]),
    )

    assert target.ewh_mm[0, 0].tolist() == [10.0, 20.0]
    assert target.month_status.tolist() == ["observed"]
