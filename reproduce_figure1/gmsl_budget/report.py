from __future__ import annotations

from pathlib import Path
from typing import Mapping, Sequence

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from .models import MonthlySeries, TrendResult


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
        "This is the standard full-ocean workflow. The budget comparison uses the fixed common 300 km coastal-buffer mask.",
        "Altimetry GIA is explicit and positive; OBD is upward-positive and is computed from the same GRACE coefficients as ocean mass.",
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
