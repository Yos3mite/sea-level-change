from __future__ import annotations

from pathlib import Path
from typing import Mapping, Sequence

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from .models import MonthlySeries, TrendResult
from .figure1_mass import prepare_figure1_series
from .trend import fit_trend


def write_diagnostics(series: Mapping[str, MonthlySeries], path: str | Path) -> None:
    figure, axis = plt.subplots(figsize=(12, 6), constrained_layout=True)
    for name, item in series.items():
        axis.plot(item.time, item.values, linewidth=1.1, label=name)
    axis.axhline(0.0, color="0.5", linewidth=0.7)
    axis.set(xlabel="Time", ylabel="Sea-level equivalent (mm)", title="GMSL budget diagnostics")
    axis.grid(alpha=0.25)
    axis.legend(ncol=2, fontsize=8)
    figure.savefig(path, dpi=180)
    plt.close(figure)


def write_mass_figure(
    series: Mapping[str, MonthlySeries],
    path: str | Path,
    hac_lags: int = 12,
    jin_target_mm_per_year: float = 1.39,
    jin_target_uncertainty_mm_per_year: float = 0.32,
) -> None:
    required = {
        "ocean_mass_csr_ssa_filled": "CSR",
        "ocean_mass_jpl_ssa_filled": "JPL",
        "ocean_mass_gsfc_ssa_filled": "GSFC",
    }
    missing = sorted(set(required) - set(series))
    if missing or "ocean_mass_ensemble" not in series:
        raise ValueError(f"mass figure is missing required series: {missing}")
    ensemble = series["ocean_mass_ensemble"]
    trend = fit_trend(ensemble, hac_lags=hac_lags)
    panel_a_members = {
        label: prepare_figure1_series(series[name], remove_trend=False, running_window=None, hac_lags=hac_lags)
        for name, label in required.items()
    }
    panel_a = prepare_figure1_series(ensemble, remove_trend=False, running_window=None, hac_lags=hac_lags)
    panel_b = prepare_figure1_series(ensemble, remove_trend=True, running_window=3, hac_lags=hac_lags)
    figure, axes = plt.subplots(2, 1, figsize=(12, 8), sharex=True, constrained_layout=True)
    colors = {"CSR": "#5B8FF9", "JPL": "#61DDAA", "GSFC": "#65789B"}
    for label, item in panel_a_members.items():
        axes[0].plot(item.time, item.values, color=colors[label], linewidth=0.9, alpha=0.65, label=label)
    axes[0].plot(panel_a.time, panel_a.values, color="#D62728", linewidth=2.0, label="Three-center mean")
    axes[0].text(
        0.015,
        0.95,
        f"Computed ensemble trend: {trend.trend_mm_per_year:.3f} ± {trend.hac_standard_error:.3f} mm/yr\n"
        f"Jin et al. reference: {jin_target_mm_per_year:.2f} ± {jin_target_uncertainty_mm_per_year:.2f} mm/yr",
        transform=axes[0].transAxes,
        va="top",
        fontsize=9,
        bbox={"facecolor": "white", "edgecolor": "0.8", "alpha": 0.9},
    )
    axes[0].set_title("Ocean-mass sea-level change: center-wise SSA, then ensemble mean")
    axes[0].set_ylabel("Deseasonalized (mm)")
    axes[0].legend(ncol=4, fontsize=8, loc="lower right")
    axes[1].plot(panel_b.time, panel_b.values, color="#D62728", linewidth=1.6)
    axes[1].axhline(0.0, color="0.5", linewidth=0.7)
    axes[1].set_title("Interannual component (linear trend removed; centered 3-month mean)")
    axes[1].set(xlabel="Time", ylabel="Anomaly (mm)")
    for axis in axes:
        axis.grid(alpha=0.25)
    figure.savefig(path, dpi=200)
    plt.close(figure)


def write_run_report(
    path: str | Path,
    run_id: str,
    config_hash: str,
    trends: Sequence[TrendResult],
    closure_available: bool,
    warnings: Sequence[str],
) -> None:
    lines = [
        f"# GMSL budget run `{run_id}`",
        "",
        f"Configuration SHA-256: `{config_hash}`",
        "",
        "## Scope",
        "",
        "The GMSL budget comparison uses a fixed common 300 km coastal-buffer mask.",
        "The Jin ocean-mass branch uses the global-ocean 1° domain without a coastal buffer; 300 km is sensitivity-only.",
        "Altimetry applies signed wet-troposphere and GIA trend corrections. Mascon product GIA removal is retained and is never applied a second time.",
        "OBD is upward-positive and is computed from the same GRACE coefficients as its spherical-harmonic mass diagnostic.",
        "",
        f"Budget closure available: **{'yes' if closure_available else 'no'}**.",
    ]
    if not closure_available:
        lines += ["", "Steric input was not supplied, so neither a steric series nor a closure residual was generated."]
    lines += ["", "## Trends", "", "| Series | Trend (mm/yr) | OLS SE | HAC SE | N |", "|---|---:|---:|---:|---:|"]
    for trend in trends:
        lines.append(
            f"| {trend.series_name} | {trend.trend_mm_per_year:.4f} | "
            f"{trend.ols_standard_error:.4f} | {trend.hac_standard_error:.4f} | {trend.n_obs} |"
        )
    lines += ["", "## Scientific warnings", ""]
    lines += [f"- {warning}" for warning in warnings] or ["- None."]
    Path(path).write_text("\n".join(lines) + "\n", encoding="utf-8")
