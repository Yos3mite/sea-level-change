# Figure 3 Regional TWS Pipeline Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build a config-driven Figure 3 pipeline that produces the paper-Mascon reproduction first and can later generate a directly comparable version from a local custom Level-3 grid product.

**Architecture:** Product adapters normalize paper Mascon and custom Level-3 files into a `MonthlyGridSeries`. A shared mask and integration layer converts each product into continental equivalent-sea-level series, while a separate temporal and plotting layer creates the four-panel Figure 3 artifact bundle. Paper reconstruction values may replace only missing observed months; every replacement and source hash is preserved.

**Tech Stack:** Python 3.12, NumPy, SciPy, netCDF4, pandas, GeoPandas, Shapely, PyProj, Matplotlib, pytest, Natural Earth Admin-0 polygons.

## Global Constraints

- Work only inside `D:\AAAA海平面变化\SaGEA-main(L2-L3)原版\SaGEA-main` for code, configuration, tests, outputs, and logs.
- Do not invoke the L2toL3 skill or read/copy code or outputs from excluded repositories.
- Preserve valid observed Mascon months; Xie–Yi values may fill only missing months through 2022-12.
- Use the canonical 1° grid with latitude `-89.5..89.5` and longitude `-179.5..179.5`.
- Use six mutually exclusive continental regions; exclude Greenland and Antarctica from continental TWS.
- Apply a land-side 300 km coastal exclusion and record the Natural Earth proxy limitation.
- Convert regional mass to equivalent global-mean sea-level millimeters using one configured global ocean area.
- Paper mode must stop on internal event-window gaps; custom mode may preserve and report NaN.
- OBD is excluded. Do not apply GIA again to products whose metadata says it is already corrected.
- Write every scientific result to a new versioned output directory; never overwrite or reinterpret an earlier result.
- Use `apply_patch` for source and documentation edits. Stage and commit only files created or modified for this feature.

---

### Task 1: Acquire and Register External Figure 3 Inputs

**Files:**
- Create: `data/external_downloads/xie_yi_2025_mascon_gapfilled/source_manifest.json`
- Create: `data/external_downloads/natural_earth/admin_0_countries/source_manifest.json`
- Create: `pysrc/reference_products/figure3/__init__.py`
- Create: `pysrc/reference_products/figure3/provenance.py`
- Test: `tests/scientific/test_figure3_provenance.py`
- Modify: `data_manifest.md`

**Interfaces:**
- Consumes: downloaded raw files from Figshare DOI `10.6084/m9.figshare.25805092.v2` and Natural Earth Admin-0 countries.
- Produces: `sha256_file(path: Path) -> str` and `build_source_manifest(source: dict, files: list[Path], root: Path) -> dict`.

- [ ] **Step 1: Download official source files without modifying them**

Use the in-app browser to open the Figshare DOI and download every CSR/JPL/GSFC gap-filled data file plus its readme into:

```text
data/external_downloads/xie_yi_2025_mascon_gapfilled/
```

Download Natural Earth Admin-0 countries at a fixed published version into:

```text
data/external_downloads/natural_earth/admin_0_countries/
```

Record the visible Figshare version, file names, byte sizes, and any publisher checksum before leaving the page. Do not rename downloaded files.

- [ ] **Step 2: Write the failing provenance test**

```python
from pathlib import Path

from pysrc.reference_products.figure3.provenance import build_source_manifest, sha256_file


def test_build_source_manifest_records_relative_paths_sizes_and_hashes(tmp_path: Path):
    raw = tmp_path / "raw"
    raw.mkdir()
    first = raw / "a.bin"
    first.write_bytes(b"abc")

    manifest = build_source_manifest(
        {"title": "source", "version": "v2", "url": "https://example.test"},
        [first],
        root=tmp_path,
    )

    assert manifest["source"]["version"] == "v2"
    assert manifest["files"] == [
        {
            "path": "raw/a.bin",
            "bytes": 3,
            "sha256": sha256_file(first),
        }
    ]
```

- [ ] **Step 3: Run the provenance test and confirm the expected failure**

Run:

```powershell
$env:PYTHONPATH="$PWD\.runtime\site-packages;$PWD"
& 'C:\Users\Alan\.cache\codex-runtimes\codex-primary-runtime\dependencies\python\python.exe' -m pytest tests/scientific/test_figure3_provenance.py -q
```

Expected: FAIL because `pysrc.reference_products.figure3.provenance` does not exist.

- [ ] **Step 4: Implement immutable-file inventory helpers**

Implement `sha256_file()` with 1 MiB streaming blocks. Implement `build_source_manifest()` so it sorts resolved files by relative POSIX path, rejects files outside `root`, and returns:

```python
{
    "source": source,
    "files": [
        {"path": relative_path, "bytes": size, "sha256": digest},
    ],
}
```

- [ ] **Step 5: Generate and validate both real source manifests**

Include DOI/URL, publication or dataset version, download timestamp, license text shown by the publisher, and file inventory. Verify every registered path exists and every SHA-256 has 64 lowercase hexadecimal characters.

- [ ] **Step 6: Update the project data manifest**

Add entries for the Xie–Yi and Natural Earth downloads. State the Xie–Yi coverage (`2002-04—2022-12`), product centers and versions, and the fact that Natural Earth is a reproducible proxy for unpublished paper polygons.

- [ ] **Step 7: Run the provenance test and commit**

Expected: PASS.

```powershell
git add -- pysrc/reference_products/figure3/__init__.py pysrc/reference_products/figure3/provenance.py tests/scientific/test_figure3_provenance.py data_manifest.md data/external_downloads/xie_yi_2025_mascon_gapfilled/source_manifest.json data/external_downloads/natural_earth/admin_0_countries/source_manifest.json
git commit -m "data: register Figure 3 source products"
```

Do not stage the large raw external files if `.gitignore` excludes them.

---

### Task 2: Define Monthly Grid Types and Product Adapters

**Files:**
- Create: `pysrc/reference_products/figure3/types.py`
- Create: `pysrc/reference_products/figure3/adapters.py`
- Test: `tests/scientific/test_figure3_adapters.py`

**Interfaces:**
- Consumes: CSR/JPL/GSFC and custom-Level-3 NetCDF paths plus per-product configuration.
- Produces: `MonthlyGridSeries`, `load_mascon(config: dict) -> MonthlyGridSeries`, and `load_custom_l3(config: dict) -> MonthlyGridSeries`.

- [ ] **Step 1: Write failing validation and adapter tests**

Use temporary NetCDF files with deliberately different dimension orders and units:

```python
def test_monthly_grid_series_rejects_duplicate_months():
    with pytest.raises(ValueError, match="duplicate months"):
        MonthlyGridSeries(
            source_id="x",
            months=np.array(["2023-01", "2023-01"]),
            lat=np.array([-0.5, 0.5]),
            lon=np.array([-0.5, 0.5]),
            ewh_mm=np.zeros((2, 2, 2)),
            valid_month=np.ones(2, dtype=bool),
            month_status=np.array(["observed", "observed"]),
            metadata={},
        )


def test_load_mascon_normalizes_cm_to_mm_and_time_lat_lon(tmp_path: Path):
    path = write_test_mascon(tmp_path, units="cm", dimension_order=("lon", "lat", "time"))
    series = load_mascon(mascon_config(path))
    assert series.ewh_mm.shape == (2, 2, 3)
    assert series.ewh_mm[0, 0, 0] == 10.0
    assert series.months.tolist() == ["2023-01", "2023-02"]
```

Also test descending latitude, `0..360` longitude conversion, fill-value-to-NaN conversion, and `valid_month` extraction for custom L3.

- [ ] **Step 2: Run the adapter tests and confirm failure**

Run the test file with the project runtime. Expected: FAIL because the types and loaders do not exist.

- [ ] **Step 3: Implement `MonthlyGridSeries` validation**

Use a frozen dataclass. Enforce unique sorted `YYYY-MM` labels, strictly increasing latitude and longitude, exact `(time, lat, lon)` field shape, finite coordinates, one status per month, and statuses limited to `observed`, `reconstructed`, or `missing`. Require missing months to contain only NaN fields.

- [ ] **Step 4: Implement configuration-driven NetCDF normalization**

The loaders must read configured variable and dimension names, decode time with `netCDF4.num2date`, transpose into `(time, lat, lon)`, convert recognized `m`, `cm`, and `mm` water-equivalent units to millimeters, normalize longitude to `[-180, 180)`, sort coordinates, and preserve correction metadata.

- [ ] **Step 5: Run tests and commit**

Expected: all adapter tests PASS.

```powershell
git add -- pysrc/reference_products/figure3/types.py pysrc/reference_products/figure3/adapters.py tests/scientific/test_figure3_adapters.py
git commit -m "feat: add Figure 3 grid product adapters"
```

---

### Task 3: Merge Observations, Reconstructions, and Target Grids

**Files:**
- Create: `pysrc/reference_products/figure3/merge.py`
- Create: `pysrc/reference_products/figure3/regrid.py`
- Test: `tests/scientific/test_figure3_merge_regrid.py`

**Interfaces:**
- Consumes: two `MonthlyGridSeries` objects for the same center.
- Produces: `align_reconstruction_baseline(observed, reconstructed, overlap_months) -> tuple[MonthlyGridSeries, dict]`, `fill_missing_months(observed, reconstructed, start, end) -> MonthlyGridSeries`, and `nearest_regrid(series, target_lat, target_lon) -> MonthlyGridSeries`.

- [ ] **Step 1: Write failing merge tests**

```python
def test_fill_missing_months_never_overwrites_observed_values():
    observed = grid_series(months=["2017-06", "2017-07"], values=[1.0, np.nan])
    reconstructed = grid_series(
        months=["2017-06", "2017-07"],
        values=[101.0, 102.0],
        status="reconstructed",
    )
    merged = fill_missing_months(observed, reconstructed, "2017-06", "2017-07")
    assert merged.ewh_mm[:, 0, 0].tolist() == [1.0, 102.0]
    assert merged.month_status.tolist() == ["observed", "reconstructed"]
```

Add tests that a reconstruction cannot extend beyond its registered end, an insufficient overlap raises, an additive bias is removed without changing differences, and nearest-neighbor regridding handles the dateline.

- [ ] **Step 2: Run the tests and confirm failure**

Expected: FAIL because merge and regrid functions do not exist.

- [ ] **Step 3: Implement additive baseline alignment**

Estimate a per-grid median `observed - reconstructed` bias over configured common observed months. Require at least 12 finite overlapping months per grid; leave cells without support as NaN. Return diagnostics containing overlap months, finite-cell fraction, bias min/mean/max, and splice RMS before and after alignment.

- [ ] **Step 4: Implement observed-first month filling**

Build the complete monthly axis, copy observed finite fields first, and fill only months whose observed status is missing. Reject any configuration that requests reconstructed months after 2022-12. Preserve per-month provenance.

- [ ] **Step 5: Implement deterministic nearest-neighbor regridding**

Use `scipy.interpolate.RegularGridInterpolator(method="nearest", bounds_error=False, fill_value=np.nan)` after coordinate normalization. Treat longitude cyclically by appending wrapped edge columns. Preserve months, validity, status, and metadata.

- [ ] **Step 6: Run tests and commit**

```powershell
git add -- pysrc/reference_products/figure3/merge.py pysrc/reference_products/figure3/regrid.py tests/scientific/test_figure3_merge_regrid.py
git commit -m "feat: merge reconstructed Figure 3 mascon grids"
```

---

### Task 4: Build Continental Masks and Regional Integration

**Files:**
- Create: `pysrc/reference_products/figure3/masks.py`
- Create: `pysrc/reference_products/figure3/integrate.py`
- Test: `tests/scientific/test_figure3_masks_integrate.py`

**Interfaces:**
- Consumes: Natural Earth shapefile, canonical grid, `MonthlyGridSeries`, and configured global ocean area.
- Produces: `ContinentMaskSet`, `build_continent_masks(...) -> ContinentMaskSet`, `cell_areas_m2(lat, lon, radius_m) -> np.ndarray`, and `integrate_regions(series, masks, ocean_area_m2) -> RegionalSeries`.

- [ ] **Step 1: Write failing analytical integration tests**

```python
def test_integrate_regions_uses_mass_not_regional_mean():
    areas = np.array([[1.0, 2.0]])
    field_mm = np.array([[[10.0, 10.0]]])
    masks = two_cell_region_masks(areas)
    result = integrate_regions(
        grid_series_from_mm(field_mm),
        masks,
        ocean_area_m2=6.0,
        water_density_kg_m3=1000.0,
    )
    assert result.values_mm["africa"][0] == pytest.approx(5.0)
```

Add tests for exact spherical cell area, mutually exclusive regions, `total == union(six regions)`, Greenland and Antarctica exclusion, and the 300 km boundary using a synthetic coastline with known point distances.

- [ ] **Step 2: Run the tests and confirm failure**

Expected: FAIL because mask and integration modules do not exist.

- [ ] **Step 3: Implement Natural Earth continent construction**

Read the configured Admin-0 shapefile with GeoPandas. Resolve both current and legacy continent field names from an explicit ordered list. Dissolve countries into six regions, remove Greenland by sovereign/admin name before North America dissolution, and discard Antarctica. Use Shapely vectorized point-in-polygon operations on canonical cell centers.

- [ ] **Step 4: Implement geodesic coastal exclusion**

Densify exterior coastline rings to at most 0.25° vertex spacing. Convert coast vertices and land cell centers to three-dimensional unit-sphere coordinates, query nearest coast vertices with `scipy.spatial.cKDTree`, convert chord distance to great-circle kilometers, and exclude land cells with distance `< 300.0 km`. Ignore interior lake rings so the rule represents ocean coastline.

- [ ] **Step 5: Implement exact spherical cell areas and ESL integration**

Use latitude-cell edges clipped at ±90°:

```python
area = radius_m**2 * delta_lon_rad * (sin(lat_north) - sin(lat_south))
regional_esl_mm = sum(ewh_mm * area) / ocean_area_m2
```

Because both EWH and ESL are in millimeters and the same water density cancels, do not introduce an unnecessary density conversion in the implementation. Store contributing area and valid-cell fraction for each region/month.

- [ ] **Step 6: Write the mask NetCDF contract**

The mask writer must include `region_id`, `land_mask`, `distance_to_coast_km`, `coastal_buffer_excluded`, and `cell_area_m2`, with region-name metadata and source hashes.

- [ ] **Step 7: Run tests and commit**

```powershell
git add -- pysrc/reference_products/figure3/masks.py pysrc/reference_products/figure3/integrate.py tests/scientific/test_figure3_masks_integrate.py
git commit -m "feat: integrate Figure 3 continental TWS"
```

---

### Task 5: Implement Temporal Processing and Paper Metrics

**Files:**
- Create: `pysrc/reference_products/figure3/temporal.py`
- Create: `pysrc/reference_products/figure3/metrics.py`
- Test: `tests/scientific/test_figure3_temporal_metrics.py`

**Interfaces:**
- Consumes: per-center `RegionalSeries` objects and event definitions.
- Produces: `process_interannual(series) -> RegionalSeries`, `combine_centers(series_by_center) -> RegionalEnsemble`, and `event_metrics(ensemble, events, paper_references) -> list[dict]`.

- [ ] **Step 1: Write failing order-of-operations tests**

Use a synthetic monthly trend plus seasonal cycle and assert that processing removes both while retaining a known interannual pulse. Test that center processing happens before averaging, centered smoothing requires three finite consecutive months, and paper mode rejects an internal event-window NaN.

- [ ] **Step 2: Write failing signed endpoint tests**

```python
def test_event_metrics_preserve_reduction_sign():
    rows = event_metrics(
        ensemble_with_total({"2014-10": 1.0, "2015-12": -5.37}),
        events=[{"id": "event", "start": "2014-10", "end": "2015-12"}],
        paper_references={"event": {"total": -6.37}},
    )
    total = next(row for row in rows if row["region"] == "total")
    assert total["change_mm"] == pytest.approx(-6.37)
    assert total["reduction_magnitude_mm"] == pytest.approx(6.37)
```

- [ ] **Step 3: Implement temporal processing**

For every center and region: subtract full-period monthly climatology, fit OLS against decimal years on finite months, subtract the fitted line, then calculate a strict 3-month centered mean. Keep centered and smoothed arrays separately.

- [ ] **Step 4: Implement center ensemble statistics**

For each month and region compute arithmetic mean, sample standard deviation when at least two centers exist, minimum, maximum, and valid-center count. Paper mode requires all configured event-window months to be finite after reconstruction and processing.

- [ ] **Step 5: Implement reference metrics without tuning**

Store start, end, signed change, reduction magnitude, paper reference, difference from paper, status, basis, and valid-center counts. Leave unreported Asia/Europe/Oceania reference values empty rather than digitizing the paper plot.

- [ ] **Step 6: Run tests and commit**

```powershell
git add -- pysrc/reference_products/figure3/temporal.py pysrc/reference_products/figure3/metrics.py tests/scientific/test_figure3_temporal_metrics.py
git commit -m "feat: process Figure 3 interannual regional series"
```

---

### Task 6: Build the Figure 3 Orchestrator and Artifact Bundle

**Files:**
- Create: `pysrc/reference_products/build_figure3_regional_tws.py`
- Test: `tests/scientific/test_figure3_pipeline.py`

**Interfaces:**
- Consumes: a complete Figure 3 JSON configuration and all modules from Tasks 1–5.
- Produces: `build_figure3(config_path: Path, project_root: Path | None = None) -> dict[str, Path]`.

- [ ] **Step 1: Write a failing end-to-end synthetic test**

Create three tiny center NetCDFs, one reconstruction file with exactly one missing-month replacement, a synthetic continent mask fixture, and a paper-mode config. Assert creation of:

```python
expected = {
    "png", "pdf", "plotting_data", "regional_by_center", "metrics",
    "masks_netcdf", "config_snapshot", "method_report", "manifest",
}
assert set(outputs) == expected
assert all(path.is_file() for path in outputs.values())
```

Assert that the manifest says one reconstructed month, no event gaps, no OBD, and no extra GIA.

- [ ] **Step 2: Run the pipeline test and confirm failure**

Expected: FAIL because the builder does not exist.

- [ ] **Step 3: Implement configuration loading and orchestration**

Validate schema sections for mode, inputs, target grid, masks, integration, time, events, references, and output. Resolve relative paths against project root. Build or load one mask set, process each center independently, combine the ensemble, calculate metrics, then plot.

- [ ] **Step 4: Implement the four-panel plot**

Use a 2×2 Matplotlib layout matching the paper structure. Use one stable color per region, a thicker Total curve, consistent colors in difference bars, zero lines, panel labels, year axes, and exact development shading. Do not apply display offsets or value scaling.

- [ ] **Step 5: Implement auditable outputs**

Write plotting CSV, per-center CSV including month provenance, metrics CSV, mask NetCDF, sorted config snapshot, method report, and manifest. The manifest must include configuration/input/output hashes, source versions, observed/reconstructed month lists, mask areas, ocean area, processing order, valid counts, warnings, and paper-reference differences.

- [ ] **Step 6: Run the synthetic pipeline test and the full new test group**

```powershell
& 'C:\Users\Alan\.cache\codex-runtimes\codex-primary-runtime\dependencies\python\python.exe' -m pytest tests/scientific/test_figure3_provenance.py tests/scientific/test_figure3_adapters.py tests/scientific/test_figure3_merge_regrid.py tests/scientific/test_figure3_masks_integrate.py tests/scientific/test_figure3_temporal_metrics.py tests/scientific/test_figure3_pipeline.py -q
```

Expected: all PASS.

- [ ] **Step 7: Commit**

```powershell
git add -- pysrc/reference_products/build_figure3_regional_tws.py tests/scientific/test_figure3_pipeline.py
git commit -m "feat: build auditable Figure 3 artifact bundle"
```

---

### Task 7: Add Paper and Custom Configurations

**Files:**
- Create: `config/figure03_paper_mascon.json`
- Create: `config/figure03_custom_l3.template.json`
- Test: `tests/scientific/test_figure3_configs.py`

**Interfaces:**
- Consumes: real downloaded and local Mascon paths.
- Produces: validated paper and custom configurations.

- [ ] **Step 1: Write a failing configuration contract test**

Assert that both configs resolve inside the project boundary except explicitly registered read-only input files, use the same grid/mask/ocean area/events, and differ only in mode/input/output sections. Assert paper mode has exactly CSR/JPL/GSFC and reconstruction inputs for the same centers.

- [ ] **Step 2: Write the paper configuration**

Register local CSR RL06.3, JPL RL06.3Mv04 CRI, GSFC RL06v2.0 and downloaded Xie–Yi files with exact variable names, dimension names, units, correction states, source versions, and hashes. Set:

```json
{
  "time": {"start": "2013-11", "end": "2024-10"},
  "grid": {"spacing_deg": 1.0},
  "mask": {"coastal_exclusion_km": 300.0},
  "integration": {"global_ocean_area_m2": 3.618e14},
  "processing": ["monthly_climatology", "ols_detrend", "centered_3_month_mean"]
}
```

Include the eight paper reference values from the design spec.

- [ ] **Step 3: Write the custom Level-3 template**

Use the existing project-default multisource Level-3 path and its concrete variable contract:

```json
{
  "input": {
    "path": "results/target_b_l2_to_l3_multisource/target_b_custom_l3_multisource_201311_202410.nc",
    "variables": {
      "time": "time",
      "lat": "lat",
      "lon": "lon",
      "field": "field",
      "land_mask": "land_mask",
      "valid_month": "valid_month"
    }
  }
}
```

Custom mode must preserve missing months and use a different output directory.

- [ ] **Step 4: Run the contract test and commit**

```powershell
git add -- config/figure03_paper_mascon.json config/figure03_custom_l3.template.json tests/scientific/test_figure3_configs.py
git commit -m "config: add paper and custom Figure 3 runs"
```

---

### Task 8: Run the Paper Reproduction and Perform Structural and Visual QA

**Files:**
- Create: `results/figure03_paper_mascon_20260820_v1/*`
- Modify: `data_manifest.md`
- Modify: `docs/target_b_l2_to_l3_method.md`

**Interfaces:**
- Consumes: `config/figure03_paper_mascon.json`.
- Produces: the final paper-Mascon Figure 3 artifact bundle and project documentation entries.

- [ ] **Step 1: Run the real paper configuration**

```powershell
$env:PYTHONPATH="$PWD\.runtime\site-packages;$PWD"
& 'C:\Users\Alan\.cache\codex-runtimes\codex-primary-runtime\dependencies\python\python.exe' -m pysrc.reference_products.build_figure3_regional_tws --config config/figure03_paper_mascon.json
```

Expected: nine outputs in `results/figure03_paper_mascon_20260820_v1/`, zero internal event-window gaps, and explicit reconstructed-month provenance.

- [ ] **Step 2: Verify numerical invariants**

Check:

```text
plotting Total equals the ensemble Total column
Total mask equals the union of six regional masks
no region masks overlap
all event endpoint values are finite
all three centers contribute to ordinary observed months
reconstruction is used only for originally missing months through 2022-12
all input and output hashes recompute exactly
```

Report signed differences from the eight paper values. Do not alter processing to reduce those differences.

- [ ] **Step 3: Re-run and verify byte-stable scientific outputs**

Run the same configuration again. Require plotting CSV, per-center CSV, metrics CSV, mask NetCDF, config snapshot, method report, and manifest scientific contents to be identical. If image/PDF metadata timestamps differ, compare rendered pixels and scientific source hashes instead of requiring byte identity.

- [ ] **Step 4: Render and inspect both image formats**

Inspect the PNG directly. Render the PDF with Poppler and inspect the rendered page. Reject clipped legends, overlapping labels, wrong panel order, missing bars, broken Chinese glyphs, inconsistent colors, or unreadable axis labels.

- [ ] **Step 5: Run the full relevant test suite**

```powershell
& 'C:\Users\Alan\.cache\codex-runtimes\codex-primary-runtime\dependencies\python\python.exe' -m pytest tests/scientific/test_figure3_*.py tests/scientific/test_figure2_custom_l3.py tests/scientific/test_figure1_custom_l3.py -q
```

Expected: zero failures.

- [ ] **Step 6: Update documentation and commit code/config/docs only**

Record the completed Figure 3 input versions, processing definition, output directory, metrics, warnings and known Natural Earth mask limitation. Do not commit ignored large source data or generated plots unless the project policy explicitly tracks them.

```powershell
git add -- data_manifest.md docs/target_b_l2_to_l3_method.md
git commit -m "docs: record paper Figure 3 reproduction"
```

---

## Completion Gate

The task is complete only when:

- the official Xie–Yi and Natural Earth inputs are downloaded and hashed;
- all new synthetic tests pass;
- the real paper configuration creates the complete nine-file artifact bundle;
- paper event windows contain no internal missing values;
- the eight signed paper-reference comparisons are reported without tuning;
- PNG and rendered PDF pass visual inspection;
- the custom-Level-3 template validates against the same regional interface;
- no excluded repository or L2toL3 workflow was accessed or invoked.
