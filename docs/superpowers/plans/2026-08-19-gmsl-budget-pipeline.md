# GMSL Budget Pipeline Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build a clean, tested, provenance-rich full-ocean GMSL processing pipeline for CMEMS/C3S altimetry, explicit GIA correction, GRACE/GRACE-FO ocean mass, and GRACE-derived ocean-bottom deformation (OBD), while leaving all historical reproduction scripts unchanged.

**Architecture:** A new `reproduce_figure1/gmsl_budget` package owns typed monthly-series models, configuration, masks, CMEMS ingestion, GIA, GRACE/GFC parsing, OBD conversion, trend fitting, provenance, and orchestration. The CLI reads one JSON configuration and writes one immutable run directory; all scientific constants and sign conventions are explicit. Production behavior is developed test-first, with small synthetic NetCDF/GFC fixtures and one final run against the local scientific data.

**Tech Stack:** Python 3.12; NumPy 1.26; Pandas 2.2; Xarray 2023.6; SciPy 1.13; Statsmodels; Matplotlib 3.8; Pytest 7.4; SaGEA local Python package; ICGEM CSR RL06/RL06.3 DDK1 spherical harmonics.

## Global Constraints

- Do not modify or import any historical scripts in `reproduce_figure1`; the new package is independent.
- Produce standard full-ocean results only; do not produce a `66°S–66°N` Jin branch.
- Use the local monthly C3S two-satellite file `D:\temp_sealevel_data\cmems_sla_monthly.nc`; do not aggregate daily altimetry.
- Keep raw and GIA-corrected altimetry series; scalar GIA defaults to `+0.30 mm/yr` and must never be applied twice.
- Use fixed spatial support. The main budget result uses the same 300 km buffered mask for altimetry, GRACE ocean mass, and OBD.
- Convert EWH to sea-level equivalent with `rho_freshwater/rho_seawater = 1000/1028`.
- Compute OBD from the same GRACE/GRACE-FO coefficients and preprocessing used for the mass term; do not hardcode `−0.11 mm/yr`.
- Treat IAP `SSL_2000m` as total steric. Do not add SIO halosteric to it. Steric remains optional in this implementation.
- Use `float64` in scientific calculations and actual calendar dates in trend time coordinates.
- Fit intercept, linear trend, annual, and semiannual terms; report OLS and HAC standard errors with default HAC lag 12.
- Preserve data source URL/DOI, version, input hashes, configuration hash, units, signs, and processing history.
- Do not commit downloaded NetCDF, GFC, ZIP, NPY, NPZ, PNG, or generated run outputs.

---

## Execution Setup: Isolated feature worktree

- [ ] **Step 1: Add and commit the worktree/data exclusions on `main`**

Create `.gitignore` with:

```gitignore
.worktrees/
__pycache__/
.pytest_cache/
*.py[cod]
*.nc
*.gfc
*.zip
*.npy
*.npz
reproduce_figure1/output/optimized_budget/
reproduce_figure1/data/external/
```

Run:

```text
git add .gitignore
git commit -m "chore: exclude generated GMSL data"
git check-ignore -q .worktrees
```

Expected: the commit contains only `.gitignore`, and `git check-ignore` exits 0.

- [ ] **Step 2: Create the isolated implementation branch and verify its baseline**

Run:

```text
git worktree add .worktrees/gmsl-budget-pipeline -b feature/gmsl-budget-pipeline
D:\python\python.exe -m pytest reproduce_figure1/tests -v
```

Expected: the worktree is created on `feature/gmsl-budget-pipeline`. If the historical project has no collected tests, Pytest exit code 5 is recorded as the empty baseline rather than described as a passing suite.

---

### Task 1: Repository hygiene, package foundation, and validated configuration

**Files:**
- Create: `reproduce_figure1/pyproject.toml`
- Create: `reproduce_figure1/gmsl_budget/__init__.py`
- Create: `reproduce_figure1/gmsl_budget/models.py`
- Create: `reproduce_figure1/gmsl_budget/config.py`
- Create: `reproduce_figure1/configs/gmsl_budget_main.json`
- Create: `reproduce_figure1/tests/gmsl_budget/test_config_models.py`

**Interfaces:**
- Produces: `MonthlySeries(time, values, name, units, metadata)`, `SpatialMask(latitude, longitude, ocean_fraction, support, metadata)`, `TrendResult`, and `PipelineConfig.load(path)`.
- Produces: `PipelineConfig.resolved_dict()` and `PipelineConfig.sha256()` for all later tasks.

- [ ] **Step 1: Write the failing configuration and model tests**

```python
def test_config_rejects_non_positive_density(tmp_path):
    path = write_config(tmp_path, rho_seawater=0.0)
    with pytest.raises(ValueError, match="rho_seawater"):
        PipelineConfig.load(path)

def test_monthly_series_rejects_duplicate_months():
    time = pd.to_datetime(["2020-01-15", "2020-01-20"])
    with pytest.raises(ValueError, match="duplicate month"):
        MonthlySeries(time, np.array([1.0, 2.0]), "x", "mm", {})

def test_config_hash_is_independent_of_json_key_order(tmp_path):
    left = PipelineConfig.load(write_config(tmp_path / "a", reverse=False))
    right = PipelineConfig.load(write_config(tmp_path / "b", reverse=True))
    assert left.sha256() == right.sha256()
```

- [ ] **Step 2: Run the tests and verify RED**

Run: `D:\python\python.exe -m pytest reproduce_figure1/tests/gmsl_budget/test_config_models.py -v`

Expected: collection fails because `gmsl_budget.models` and `gmsl_budget.config` do not exist.

- [ ] **Step 3: Implement immutable dataclasses and JSON validation**

```python
@dataclass(frozen=True)
class MonthlySeries:
    time: pd.DatetimeIndex
    values: np.ndarray
    name: str
    units: str
    metadata: Mapping[str, Any]

    def __post_init__(self) -> None:
        months = self.time.to_period("M")
        if months.duplicated().any():
            raise ValueError("duplicate month in MonthlySeries")
        if np.asarray(self.values, dtype=np.float64).shape != (len(self.time),):
            raise ValueError("values must have shape (time,)")
```

`PipelineConfig.load()` must validate required paths, positive densities, `0 < min_valid_weight_fraction <= 1`, `gia_mode in {"scalar", "spatial_caron"}`, `hac_lags >= 0`, and OBD/mass preprocessing equality.

- [ ] **Step 4: Run the focused tests and then the package test directory**

Run: `D:\python\python.exe -m pytest reproduce_figure1/tests/gmsl_budget/test_config_models.py -v`

Expected: all tests pass.

- [ ] **Step 5: Commit the task**

```bash
git add reproduce_figure1/pyproject.toml reproduce_figure1/gmsl_budget reproduce_figure1/configs reproduce_figure1/tests/gmsl_budget/test_config_models.py
git commit -m "feat: add validated GMSL pipeline configuration"
```

### Task 2: Calendar-aware trend and uncertainty estimation

**Files:**
- Create: `reproduce_figure1/gmsl_budget/trend.py`
- Create: `reproduce_figure1/tests/gmsl_budget/test_trend.py`

**Interfaces:**
- Consumes: `MonthlySeries` and `TrendResult` from Task 1.
- Produces: `fit_trend(series: MonthlySeries, hac_lags: int = 12) -> TrendResult`.

- [ ] **Step 1: Write failing synthetic trend tests**

```python
def test_fit_trend_recovers_calendar_trend_and_harmonics():
    time = pd.date_range("2010-01-15", periods=120, freq="MS") + pd.Timedelta(days=14)
    year = decimal_year(time)
    values = 7.0 + 3.2 * (year - year.mean()) + 4.0 * np.sin(2*np.pi*year) - 1.5 * np.cos(4*np.pi*year)
    result = fit_trend(MonthlySeries(time, values, "synthetic", "mm", {}), 12)
    assert result.trend_mm_per_year == pytest.approx(3.2, abs=1e-10)
    assert result.n_obs == 120

@pytest.mark.parametrize("periods", [24, 35])
def test_fit_trend_rejects_short_series(periods):
    series = MonthlySeries(pd.date_range("2020-01-15", periods=periods, freq="MS"), np.ones(periods), "x", "mm", {})
    with pytest.raises(ValueError, match="36 valid months"):
        fit_trend(series)
```

- [ ] **Step 2: Run the tests and verify RED**

Run: `D:\python\python.exe -m pytest reproduce_figure1/tests/gmsl_budget/test_trend.py -v`

Expected: import failure for `gmsl_budget.trend`.

- [ ] **Step 3: Implement decimal-year design matrix and OLS/HAC results**

```python
def fit_trend(series: MonthlySeries, hac_lags: int = 12) -> TrendResult:
    ok = np.isfinite(series.values)
    time = series.time[ok]
    if ok.sum() < 36 or decimal_year(time).ptp() < 3.0:
        raise ValueError("trend requires at least 36 valid months spanning 3 years")
    t = decimal_year(time)
    centered = t - t.mean()
    x = np.column_stack((np.ones(len(t)), centered,
                         np.sin(2*np.pi*t), np.cos(2*np.pi*t),
                         np.sin(4*np.pi*t), np.cos(4*np.pi*t)))
    model = statsmodels.api.OLS(series.values[ok], x).fit()
    robust = model.get_robustcov_results(cov_type="HAC", maxlags=hac_lags)
```

Reject rank-deficient or ill-conditioned design matrices. Store trend, OLS/HAC standard errors, residual lag-1 correlation, number of observations, date range, and missing months.

- [ ] **Step 4: Run RED→GREEN verification and the full current suite**

Run: `D:\python\python.exe -m pytest reproduce_figure1/tests/gmsl_budget -v`

Expected: all current tests pass without warnings.

- [ ] **Step 5: Commit the task**

```bash
git add reproduce_figure1/gmsl_budget/trend.py reproduce_figure1/tests/gmsl_budget/test_trend.py
git commit -m "feat: add calendar-aware harmonic trend fitting"
```

### Task 3: Fixed masks, coastal fractions, and CMEMS monthly GMSL

**Files:**
- Create: `reproduce_figure1/gmsl_budget/masks.py`
- Create: `reproduce_figure1/gmsl_budget/cmems.py`
- Create: `reproduce_figure1/tests/gmsl_budget/test_masks_cmems.py`
- Create: `reproduce_figure1/tests/fixtures/make_cmems_fixture.py`

**Interfaces:**
- Consumes: `SpatialMask`, `MonthlySeries`, and configuration values.
- Produces: `cell_area_weights(latitude, longitude) -> np.ndarray`.
- Produces: `load_cdt_ocean_fraction(land_mask_mat, target_lat, target_lon) -> xr.DataArray`.
- Produces: `load_cdt_coast_distance(distance_mat, target_lat, target_lon) -> xr.DataArray`.
- Produces: `build_altimetry_mask(sla: xr.DataArray, ocean_fraction: xr.DataArray) -> SpatialMask`.
- Produces: `buffer_ocean_mask(mask: SpatialMask, coast_distance_km: xr.DataArray, distance_km: float) -> SpatialMask`.
- Produces: `read_cmems_gmsl(path, mask, min_valid_weight_fraction) -> MonthlySeries`.

- [ ] **Step 1: Write failing mask and area-average tests**

```python
def test_area_weights_use_cell_bounds_not_only_cosine_centers():
    lat = np.array([-60.0, 0.0, 60.0])
    lon = np.array([0.0, 120.0, 240.0])
    weights = cell_area_weights(lat, lon)
    assert weights[1, 0] > weights[0, 0]
    assert weights[:, 0].sum() == pytest.approx(4*np.pi*EARTH_RADIUS_M**2 / 3, rel=1e-12)

def test_monthly_missing_cell_does_not_renormalize_fixed_support(tmp_path):
    path, mask = make_two_month_fixture(tmp_path, second_month_has_missing_cell=True)
    result = read_cmems_gmsl(path, mask, min_valid_weight_fraction=0.995)
    assert np.isfinite(result.values[0])
    assert np.isnan(result.values[1])

def test_buffer_keeps_only_ocean_cells_at_least_300_km_from_coast():
    mask = make_three_cell_ocean_mask()
    distance = xr.DataArray([[150.0, 300.0, 350.0]], dims=("latitude", "longitude"))
    buffered = buffer_ocean_mask(mask, distance, distance_km=300.0)
    assert buffered.support.tolist() == [[False, True, True]]

def test_cdt_subcell_aggregation_returns_fractional_coastal_cell(tmp_path):
    land_path = write_cdt_fixture(tmp_path, land=np.array([[0, 1], [0, 1]]))
    fraction = load_cdt_ocean_fraction(land_path, np.array([0.0]), np.array([0.0]))
    assert fraction.item() == pytest.approx(0.5)
```

- [ ] **Step 2: Run the tests and verify RED**

Run: `D:\python\python.exe -m pytest reproduce_figure1/tests/gmsl_budget/test_masks_cmems.py -v`

Expected: imports fail because `masks.py` and `cmems.py` do not exist.

- [ ] **Step 3: Implement exact spherical cell areas and fixed-support averaging**

Use latitude/longitude bounds inferred from monotonic centers. Compute cell area as `R² × Δλ × (sin φ_north − sin φ_south)`. The monthly numerator uses only values on the fixed support, and the denominator is always the complete fixed-support weight sum; a month below the valid-weight threshold becomes NaN and records a warning.

- [ ] **Step 4: Implement CDT 1/8° ocean fractions and 300 km buffering**

Load `主要资料整理\风应力计算\wind\cdt\cdt_data\land_mask.mat`, whose `land` grid is a 1/8° binary mask, and aggregate the four 1/8° subcell samples inside each 0.25° C3S cell to obtain fractions in `{0, 0.25, 0.5, 0.75, 1}`. Load `distance2coast.mat`, whose `D` variable is documented as Haversine great-circle distance in kilometres, and sample it at target cell centres. Keep ocean cells with `D >= 300 km`, preserve the fractional ocean weights, and record the CDT citation (Greene et al., 2019, DOI `10.1029/2019GC008392`) and source hashes.

- [ ] **Step 5: Run focused and full tests**

Run: `D:\python\python.exe -m pytest reproduce_figure1/tests/gmsl_budget -v`

Expected: all tests pass.

- [ ] **Step 6: Commit the task**

```bash
git add reproduce_figure1/gmsl_budget/masks.py reproduce_figure1/gmsl_budget/cmems.py reproduce_figure1/tests
git commit -m "feat: add fixed-mask CMEMS GMSL processing"
```

### Task 4: Explicit and non-repeatable GIA correction

**Files:**
- Create: `reproduce_figure1/gmsl_budget/gia.py`
- Create: `reproduce_figure1/tests/gmsl_budget/test_gia.py`

**Interfaces:**
- Consumes: `MonthlySeries`.
- Produces: `apply_scalar_gia(series, rate_mm_per_year=0.30, reference_time=None) -> MonthlySeries`.
- Produces: `area_average_spatial_gia(gia_rate_grid, mask) -> float`.

- [ ] **Step 1: Write failing sign and duplicate-application tests**

```python
def test_positive_gia_rate_adds_three_mm_over_ten_years():
    time = pd.to_datetime(["2010-01-15", "2020-01-15"])
    raw = MonthlySeries(time, np.zeros(2), "gmsl_raw", "mm", {"gia_corrected": False})
    corrected = apply_scalar_gia(raw, 0.30, reference_time=time[0])
    assert corrected.values[1] - corrected.values[0] == pytest.approx(3.0, rel=2e-4)

def test_gia_cannot_be_applied_twice():
    corrected = make_series(metadata={"gia_corrected": True})
    with pytest.raises(ValueError, match="already GIA corrected"):
        apply_scalar_gia(corrected, 0.30)
```

- [ ] **Step 2: Run tests and verify RED**

Run: `D:\python\python.exe -m pytest reproduce_figure1/tests/gmsl_budget/test_gia.py -v`

Expected: import failure for `gmsl_budget.gia`.

- [ ] **Step 3: Implement scalar and spatial-average GIA modes**

The scalar correction is `rate × (decimal_year(time) − decimal_year(reference_time))`. Set metadata keys `gia_corrected=True`, `gia_mode`, `gia_rate_mm_per_year`, `gia_reference_time`, and `parent_series`. Reject official CMEMS indicator variables whose attributes/name say `GIA_corrected`.

- [ ] **Step 4: Run all tests and commit**

Run: `D:\python\python.exe -m pytest reproduce_figure1/tests/gmsl_budget -v`

Expected: all tests pass.

```bash
git add reproduce_figure1/gmsl_budget/gia.py reproduce_figure1/tests/gmsl_budget/test_gia.py
git commit -m "feat: make GIA correction explicit and sign-safe"
```

### Task 5: GRACE monthly adapters, month alignment, and density conversion

**Files:**
- Create: `reproduce_figure1/gmsl_budget/grace.py`
- Create: `reproduce_figure1/tests/gmsl_budget/test_grace.py`

**Interfaces:**
- Consumes: `MonthlySeries` and `SpatialMask`.
- Produces: `ewh_to_sea_level(ewh, rho_freshwater=1000.0, rho_seawater=1028.0)`.
- Produces: `align_common_months(series: Sequence[MonthlySeries]) -> list[MonthlySeries]`.
- Produces: `ensemble_mean(series_by_center: Mapping[str, MonthlySeries]) -> MonthlySeries`.
- Produces: `read_mascon_ocean_series(path, center, mask, variable=None) -> MonthlySeries`.

- [ ] **Step 1: Write failing behavior tests**

```python
def test_ewh_density_conversion_uses_fresh_to_seawater_ratio():
    assert ewh_to_sea_level(np.array([1028.0]))[0] == pytest.approx(1000.0)

def test_alignment_uses_year_month_keys_not_array_positions():
    a = make_series(["2020-01-15", "2020-03-15"], [1.0, 3.0], "CSR")
    b = make_series(["2020-02-15", "2020-03-16"], [2.0, 4.0], "JPL")
    aa, bb = align_common_months([a, b])
    assert list(aa.time.to_period("M").astype(str)) == ["2020-03"]
    assert aa.values.tolist() == [3.0]
    assert bb.values.tolist() == [4.0]

def test_ensemble_metadata_preserves_actual_centers():
    result = ensemble_mean({"CSR": csr, "GFZ": gfz, "JPL": jpl})
    assert result.metadata["centers"] == ["CSR", "GFZ", "JPL"]
```

- [ ] **Step 2: Run and verify RED**

Run: `D:\python\python.exe -m pytest reproduce_figure1/tests/gmsl_budget/test_grace.py -v`

Expected: import failure for `gmsl_budget.grace`.

- [ ] **Step 3: Implement adapters and explicit metadata checks**

Never interpolate missing months. Require center, solution/release, GIA model, filter, leakage method, density ratio, mask hash, and units in metadata. Reject a center label that disagrees with NetCDF global attributes or filename tokens.

- [ ] **Step 4: Run all tests and commit**

Run: `D:\python\python.exe -m pytest reproduce_figure1/tests/gmsl_budget -v`

Expected: all tests pass.

```bash
git add reproduce_figure1/gmsl_budget/grace.py reproduce_figure1/tests/gmsl_budget/test_grace.py
git commit -m "feat: add traceable GRACE monthly adapters"
```

### Task 6: Secure ICGEM acquisition and GFC parsing

**Files:**
- Create: `reproduce_figure1/gmsl_budget/icgem.py`
- Create: `reproduce_figure1/scripts/download_icgem_csr.py`
- Create: `reproduce_figure1/tests/gmsl_budget/test_icgem.py`
- Create: `reproduce_figure1/tests/fixtures/csr_low_degree_sample.gfc`

**Interfaces:**
- Produces: `parse_gfc(path: Path, lmax: int = 60) -> GfcEpoch`.
- Produces: `discover_zip_links(series_pages: Sequence[str], filter_name="DDK1") -> list[DownloadSpec]`.
- Produces: `download_and_extract(spec, destination, expected_sha256=None) -> list[Path]`.
- The script accepts `--destination`, `--start 2013-11`, `--end 2024-10`, and `--dry-run`.

- [ ] **Step 1: Write failing parser and download-safety tests**

```python
def test_parse_gfc_reads_dates_normalization_and_coefficients(fixture_path):
    epoch = parse_gfc(fixture_path, lmax=2)
    assert epoch.start == date(2020, 1, 1)
    assert epoch.end == date(2020, 1, 31)
    assert epoch.c[0, 0] == 0.0
    assert epoch.c[2, 0] == pytest.approx(-4.84e-4)
    assert epoch.normalization == "fully_normalized"

def test_download_rejects_html_saved_as_gfc(tmp_path):
    response = FakeResponse(b"<html>error</html>", content_type="text/html")
    with pytest.raises(ValueError, match="not a GFC or ZIP"):
        download_and_extract(make_spec(response), tmp_path)
```

- [ ] **Step 2: Run tests and verify RED**

Run: `D:\python\python.exe -m pytest reproduce_figure1/tests/gmsl_budget/test_icgem.py -v`

Expected: import failure for `gmsl_budget.icgem`.

- [ ] **Step 3: Implement secure discovery/download and deterministic parsing**

Use HTTPS only and the system CA store. Do not disable certificate or hostname verification. Download to `*.part`, validate HTTP content type/ZIP magic/GFC header, compute SHA-256, then atomically rename. Select CSR GRACE RL06 and CSR GRACE-FO RL06.3 DDK1 files whose epoch midpoint falls within 2013-11 through 2024-10. Record ICGEM page URL, DOI, license, retrieval time, and hashes in `download_manifest.json`.

- [ ] **Step 4: Verify the CLI dry run against official ICGEM pages**

Run: `D:\python\python.exe reproduce_figure1/scripts/download_icgem_csr.py --destination D:\temp_sealevel_data\grace_csr_ddk1 --start 2013-11 --end 2024-10 --dry-run`

Expected: lists GRACE and GRACE-FO archives/files, performs no writes, and exits 0.

- [ ] **Step 5: Run all tests and commit**

Run: `D:\python\python.exe -m pytest reproduce_figure1/tests/gmsl_budget -v`

Expected: all tests pass.

```bash
git add reproduce_figure1/gmsl_budget/icgem.py reproduce_figure1/scripts/download_icgem_csr.py reproduce_figure1/tests
git commit -m "feat: add secure ICGEM CSR acquisition"
```

### Task 7: Shared GRACE preprocessing and SaGEA OBD conversion

**Files:**
- Create: `reproduce_figure1/gmsl_budget/obd.py`
- Create: `reproduce_figure1/tests/gmsl_budget/test_obd.py`
- Create: `reproduce_figure1/tests/gmsl_budget/test_sagea_obd_integration.py`

**Interfaces:**
- Consumes: parsed `GfcEpoch` objects, low-degree files, Caron GIA, `SpatialMask`, and SaGEA path.
- Produces: `preprocess_grace_coefficients(epochs, config) -> ProcessedCoefficients` with a preprocessing hash.
- Produces: `convert_coefficients_to_ewh(processed) -> xr.DataArray`.
- Produces: `convert_coefficients_to_vertical_displacement(processed) -> xr.DataArray`.
- Produces: `area_average_obd(displacement, mask, preprocessing_hash) -> MonthlySeries`.

- [ ] **Step 1: Write failing Love-number scaling, sign, and shared-hash tests**

```python
def test_vertical_displacement_uses_h_over_one_plus_k():
    c = single_degree_coefficient(degree=2, value=1e-10)
    got = vertical_conversion_weights(lmax=2, k=np.array([0.0, 0.0, -0.3]), h=np.array([0.0, 0.0, -0.6]))
    assert got[2] == pytest.approx(EARTH_RADIUS_M * (-0.6) / 0.7)

def test_downward_ocean_bottom_motion_is_negative_obd():
    displacement = xr.DataArray([[[-0.002, -0.002]]], dims=("time", "latitude", "longitude"))
    result = area_average_obd(displacement, two_cell_ocean_mask(), "abc123")
    assert result.values.tolist() == pytest.approx([-2.0])

def test_mass_and_obd_reject_different_preprocessing_hashes():
    with pytest.raises(ValueError, match="preprocessing hash"):
        assert_shared_preprocessing(mass_hash="aaa", obd_hash="bbb")
```

- [ ] **Step 2: Run tests and verify RED**

Run: `D:\python\python.exe -m pytest reproduce_figure1/tests/gmsl_budget/test_obd.py -v`

Expected: import failure for `gmsl_budget.obd`.

- [ ] **Step 3: Implement one coefficient preprocessing function**

Apply low-degree replacement (degree 1, C20, C30), subtract the Caron 2018 GIA trend at each epoch midpoint, remove the temporal mean, and retain the ICGEM DDK1 filter without applying a second Gaussian filter. Serialize every setting and input hash into one preprocessing hash consumed by both mass and OBD branches.

- [ ] **Step 4: Implement EWH and vertical-displacement branches**

Use SaGEA `ConvertSHC` with `PhysicalDimensions.EWH` for the mass grid and `PhysicalDimensions.VerticalDisplacement` for OBD. Load Wang `k_l` and `h_l` from the same `LoveNumber.mat`. Synthesize both fields on the configured grid and apply identical `budget_common` weights. Keep vertical displacement in metres until the final conversion to millimetres.

- [ ] **Step 5: Run numerical integration tests**

Run: `D:\python\python.exe -m pytest reproduce_figure1/tests/gmsl_budget/test_sagea_obd_integration.py -v`

Expected: zero load returns exactly zero; a synthetic degree-2 coefficient matches the independently calculated Love-number scale; all arrays are float64 and finite.

- [ ] **Step 6: Run all tests and commit**

Run: `D:\python\python.exe -m pytest reproduce_figure1/tests/gmsl_budget -v`

Expected: all tests pass.

```bash
git add reproduce_figure1/gmsl_budget/obd.py reproduce_figure1/tests/gmsl_budget/test_obd.py reproduce_figure1/tests/gmsl_budget/test_sagea_obd_integration.py
git commit -m "feat: compute OBD from shared GRACE coefficients"
```

### Task 8: Provenance, immutable outputs, pipeline orchestration, and CLI

**Files:**
- Create: `reproduce_figure1/gmsl_budget/provenance.py`
- Create: `reproduce_figure1/gmsl_budget/pipeline.py`
- Create: `reproduce_figure1/gmsl_budget/report.py`
- Create: `reproduce_figure1/run_gmsl_budget.py`
- Create: `reproduce_figure1/tests/gmsl_budget/test_pipeline.py`

**Interfaces:**
- Produces: `sha256_file(path) -> str`.
- Produces: `run_pipeline(config: PipelineConfig) -> Path`.
- Produces: `write_run_outputs(run, output_dir) -> None`.
- CLI: `D:\python\python.exe reproduce_figure1/run_gmsl_budget.py --config reproduce_figure1/configs/gmsl_budget_main.json`.

- [ ] **Step 1: Write failing end-to-end fixture tests**

```python
def test_pipeline_without_steric_writes_valid_partial_budget(tmp_path):
    config = fixture_pipeline_config(tmp_path, steric_path=None)
    run_dir = run_pipeline(config)
    ds = xr.open_dataset(run_dir / "monthly_budget.nc")
    assert "gmsl_raw" in ds
    assert "gmsl_gia_corrected" in ds
    assert "ocean_mass" in ds
    assert "obd" in ds
    assert "steric" not in ds
    assert json.loads((run_dir / "provenance.json").read_text())["budget_closure_available"] is False

def test_pipeline_refuses_to_overwrite_different_config(tmp_path):
    first = fixture_pipeline_config(tmp_path, gia_rate=0.30)
    second = fixture_pipeline_config(tmp_path, gia_rate=0.25, run_id=first.run_id)
    run_pipeline(first)
    with pytest.raises(FileExistsError, match="different configuration hash"):
        run_pipeline(second)
```

- [ ] **Step 2: Run tests and verify RED**

Run: `D:\python\python.exe -m pytest reproduce_figure1/tests/gmsl_budget/test_pipeline.py -v`

Expected: imports fail because pipeline modules do not exist.

- [ ] **Step 3: Implement orchestration and immutable run output**

Write `monthly_budget.nc`, `monthly_budget.csv`, `trend_summary.csv`, both mask NetCDF files, `provenance.json`, `config_resolved.json`, `diagnostics.png`, and `run_report.md`. NetCDF variables must carry units, signs, source, GIA status, mask hash, and preprocessing hash. Omit steric and closure variables when no qualified steric series exists.

- [ ] **Step 4: Run fixture pipeline twice for reproducibility**

Run: `D:\python\python.exe -m pytest reproduce_figure1/tests/gmsl_budget/test_pipeline.py -v`

Expected: both runs produce identical numerical arrays and configuration hashes; all tests pass.

- [ ] **Step 5: Run the full test suite and commit**

Run: `D:\python\python.exe -m pytest reproduce_figure1/tests/gmsl_budget -v`

Expected: all tests pass with zero failures.

```bash
git add reproduce_figure1/gmsl_budget reproduce_figure1/run_gmsl_budget.py reproduce_figure1/tests/gmsl_budget/test_pipeline.py
git commit -m "feat: orchestrate auditable GMSL budget runs"
```

### Task 9: Scientific data run, external OBD validation, and final documentation

**Files:**
- Create: `reproduce_figure1/scripts/download_obd_validation.py`
- Create: `reproduce_figure1/tests/gmsl_budget/test_obd_validation.py`
- Create: `reproduce_figure1/GMSL_BUDGET_README.md`
- Modify: `reproduce_figure1/configs/gmsl_budget_main.json`

**Interfaces:**
- Produces: `compare_obd(reference: MonthlySeries, computed: MonthlySeries) -> ValidationMetrics`.
- Produces final run directory under `reproduce_figure1/output/optimized_budget/<run_id>/`.

- [ ] **Step 1: Write failing reference-alignment and metric tests**

```python
def test_obd_validation_aligns_months_and_removes_only_common_mean():
    reference = make_series(["2020-01", "2020-03"], [1.0, 3.0], "reference")
    computed = make_series(["2020-02", "2020-03"], [8.0, 10.0], "computed")
    metrics = compare_obd(reference, computed)
    assert metrics.n_common == 1
    assert metrics.common_months == ["2020-03"]
```

- [ ] **Step 2: Run the test and verify RED**

Run: `D:\python\python.exe -m pytest reproduce_figure1/tests/gmsl_budget/test_obd_validation.py -v`

Expected: `compare_obd` is missing.

- [ ] **Step 3: Implement and run the optional external validator**

The download script stores the Adhikari/Harvard Dataverse DOI, resolved file URL, license, retrieval date, and SHA-256. Validation reports common dates, mean offset, trend difference, correlation, and RMS difference; it never replaces the computed OBD series.

- [ ] **Step 4: Download official CSR DDK1 data and run the production pipeline**

Run:

```text
D:\python\python.exe reproduce_figure1/scripts/download_icgem_csr.py --destination D:\temp_sealevel_data\grace_csr_ddk1 --start 2013-11 --end 2024-10
D:\python\python.exe reproduce_figure1/run_gmsl_budget.py --config reproduce_figure1/configs/gmsl_budget_main.json
```

Expected: one immutable run directory containing all outputs listed in the design; no historical file changes.

- [ ] **Step 5: Execute scientific acceptance checks**

Run: `D:\python\python.exe -m pytest reproduce_figure1/tests/gmsl_budget -v`

Then inspect `trend_summary.csv` and require:

- C3S raw full-ocean trend is in the diagnostic range 4.08–4.13 mm/yr, or the report explains a mask/version difference;
- GIA-corrected trend minus raw trend is 0.30 mm/yr within numerical tolerance;
- OBD trend outside −0.30 to +0.10 mm/yr raises a scientific warning without changing the value;
- official CMEMS indicator comparison records that the indicator is GIA corrected, seasonally adjusted, and six-month filtered;
- no closure result exists when steric input is absent.

- [ ] **Step 6: Document commands, scientific conventions, limitations, and private-repository authentication**

`GMSL_BUDGET_README.md` must state the input product identifiers, full-ocean-only scope, GIA sign, OBD sign, fixed masks, density ratio, trend equation, output schema, and exact rerun commands. It must explain that GitHub push requires either Git Credential Manager login or a fine-grained PAT with `Contents: Read and write` for `Yos3mite/sea-level-change`.

- [ ] **Step 7: Run final verification and commit**

Run:

```text
D:\python\python.exe -m pytest reproduce_figure1/tests/gmsl_budget -v
D:\python\python.exe -m compileall -q reproduce_figure1/gmsl_budget reproduce_figure1/run_gmsl_budget.py reproduce_figure1/scripts
git status --short
```

Expected: zero test failures, compile exit 0, and only intended source/documentation changes before commit.

```bash
git add reproduce_figure1/GMSL_BUDGET_README.md reproduce_figure1/configs/gmsl_budget_main.json reproduce_figure1/scripts/download_obd_validation.py reproduce_figure1/tests/gmsl_budget/test_obd_validation.py
git commit -m "docs: validate and document optimized GMSL workflow"
```

### Task 10: Branch completion and private GitHub handoff

**Files:**
- No production file changes.

**Interfaces:**
- Consumes: completed feature branch and fresh verification output.
- Produces: a verified local branch ready for authenticated push.

- [ ] **Step 1: Re-run complete verification from the feature worktree**

Run:

```text
D:\python\python.exe -m pytest reproduce_figure1/tests/gmsl_budget -v
D:\python\python.exe -m compileall -q reproduce_figure1/gmsl_budget reproduce_figure1/run_gmsl_budget.py reproduce_figure1/scripts
git diff --check main...HEAD
git log --oneline --decorate main..HEAD
```

Expected: tests and compilation pass, `git diff --check` produces no output, and the log contains the task commits.

- [ ] **Step 2: Verify private-repository authentication without prompting**

Run: `git -c credential.interactive=never ls-remote origin HEAD`

Expected after authentication: exit 0 and a HEAD reference. Before authentication, report the credential blocker and do not claim that remote upload succeeded.

- [ ] **Step 3: Push only after authentication and remote-history reconciliation**

If the remote is empty, run `git push -u origin feature/gmsl-budget-pipeline`. If it has existing commits, fetch first and integrate without force-push. Never use `--force` or overwrite remote history.
