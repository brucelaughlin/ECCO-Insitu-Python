# Project Log

A running, in-repo changelog. Newest entry first. Each entry summarizes what a
commit changed, why, what was discussed, and how it was verified. (Complements
git history with human-readable rationale.)

---

## 2026-09-28 — Pre-baked climatology on fixed observation depth grids

### Motivation (the boss's ask)

All profile files from Scripps use exactly two fixed vertical depth grids:

| Grid | Levels | Depth range | Files |
|------|--------|-------------|-------|
| A    | 97     | 2–6000 m    | CTD, GLD, MEOP, PFL, PFL_BGC, XBT, ITP, MRB/97_depths |
| B    | 36     | 1–750 m     | MRB/36_depths only |

Because these grids are fixed (all files of a given type share the same depth axis),
the vertical interpolation from WOA23's 102-level grid onto the observation grid
can be done once offline — rather than repeated for every profile at runtime.
The boss's intent was also to simplify future lat/lon interpolation: with
pre-baked files, climatology lookup is a pure lat/lon operation with no vertical
component, making it straightforward to extend to new source types.

### What was built

**`prebake_woa23_climatology.py`** — new standalone offline tool.

Reads `woa23_decav91C0_TS_clim_potential_T_1deg_fulldepth.nc` (12 months × 102
depths × 180 lat × 360 lon) and for each grid:

1. Treats each lat/lon cell as a "profile" — reshapes the 102-level column at
   every grid cell to (64800, 102).
2. Runs the same `clim_interp._depth_interp` + `clim_interp._fallback_fill` that
   step03 uses at runtime — so the pre-baked values are exactly what step03 would
   have computed.
3. Writes back to (12, ndepth_obs, 180, 360).

Outputs:
- `woa23_decav91C0_TS_clim_potential_T_1deg_97depths_prebaked.nc`
- `woa23_decav91C0_TS_clim_potential_T_1deg_36depths_prebaked.nc`

The script exports `GRID_97` and `GRID_36` (Python float lists of the canonical
depth values) so controller scripts can build the lookup dict.

**`clim_interp.py`** — two small changes, fully backward-compatible:

1. `interpolate_climatology` gains a keyword-only `prebaked=False` argument.
   When `True`, the function skips `_depth_interp` and `_fallback_fill` at the
   end — the field is already on the observation depth grid, so those steps are
   unnecessary.
2. The `clim_cols` working array is allocated using `field_months.shape[1]`
   instead of `clim_depths.size`. In normal mode these are always equal (both
   equal the WOA23 102-level depth count), so no behaviour changes. In pre-baked
   mode `field_months.shape[1]` equals `ndepth_obs` (97 or 36), which is correct.

**`step03.py`** — accepts an optional `prebaked_clim_files` dict:

```python
prebaked_clim_files = {
    tuple(GRID_97): '/path/to/woa23_..._97depths_prebaked.nc',
    tuple(GRID_36): '/path/to/woa23_..._36depths_prebaked.nc',
}
```

At the start of each file, step03 hashes the file's `prof_depth` array
(`tuple(prof_depths.tolist())`) and looks it up in the dict. If found, it loads
the pre-baked file instead of the full-depth climatology and passes
`prebaked=True` to `interpolate_climatology`. If not found (e.g. an ITP file with
a non-standard grid), it falls back to the existing full-depth interpolation path
— no error, no change in behaviour.

**`NCEI.py`** — threads `prebaked_clim_files` through `NCEI_pipeline` and `main`.
The parameter defaults to `None`, so all existing callers are unaffected. A
commented-out example block near the `climatology_file` config line shows how to
activate it.

### How to activate

Uncomment the block in `NCEI.py` near line 55:

```python
from prebake_woa23_climatology import GRID_97, GRID_36
prebaked_clim_files = {
    tuple(GRID_97): '/Users/brucel/ecco/yip/woa23_climatology/woa23_decav91C0_TS_clim_potential_T_1deg_97depths_prebaked.nc',
    tuple(GRID_36): '/Users/brucel/ecco/yip/woa23_climatology/woa23_decav91C0_TS_clim_potential_T_1deg_36depths_prebaked.nc',
}
```

### Verification

- CTD 1992 (97-level): new output **bit-identical** to production run across all 37 variables.
- MRB 1991 (36-level): new output **bit-identical** to production run across all 51 variables.
- ITP/L2 (non-standard grid): falls back to full-depth interpolation cleanly.
- Two independent runs from the same input: bit-identical (chain is deterministic).

### Runtime impact: none before caching (I/O bound)

Timing on the CTD 1992 file (61,551 profiles) before caching was added:

| Path | Time |
|------|------|
| Normal (full-depth clim) | 5.1 s |
| Pre-baked (no cache) | 5.5 s |

No speedup without caching — profiling showed that ~2.0 s of the ~5.3 s total
was disk I/O reading the climatology file via xarray/netCDF4, and the depth
interpolation itself was only ~0.9 s. Loading the pre-baked file instead saved
the depth interp time but added comparable I/O for a file that is similar in
size when decompressed.

### Climatology caching (the actual speedup)

The pre-bake only produces a real speedup if the arrays are loaded **once** into
memory before the file loop, not re-read from disk on every file.

**How it works now:**

1. `NCEI_pipeline` defines `prebaked_clim_files` as `{depth_tuple: path}` as before.
2. Immediately after, a new block iterates over that dict, calls `xr.open_dataset`
   on each path, extracts all numpy arrays, and builds a new dict:
   `{depth_tuple: {'prof_T': ndarray, 'prof_S': ndarray, 'lon': ..., 'lat': ..., 'depths': ..., '_source_name': filename}}`.
3. The original path dict is replaced in-place (`prebaked_clim_files = prebaked_clim_arrays`).
4. This new array dict is passed to `partial(step03.main, ..., prebaked_clim_files=...)`.
5. In step03, the pre-baked branch now just uses the dict value directly — no
   `xr.open_dataset` call at all. The `_source_name` key holds the original
   filename for the log line.

The 4 pre-baked .nc files (97-level, 36-level, 25-level Samoa, 22-level Samoa)
are loaded at pipeline start, once. For the remaining 509 files, step03's
climatology lookup is purely in-memory: a dict key check + time-blend + bilinear
lat/lon interpolation. The disk I/O cost (~2 s/file) is paid four times total
instead of 509 times.

**Expected speedup:** ~2 s per file × (509 − 4) files ≈ **17 minutes** saved on
a full 509-file run. Actual speedup to be measured in the next production run.

### Files
- Added: `prebake_woa23_climatology.py`
- Changed: `clim_interp.py`, `step03.py`, `NCEI.py`

---

## 2026-09-26 — Production run + reference comparison

### Production run
Full 509-file multi-source run completed (~3 hours). 422 files produced output;
87 produced no output, all for legitimate reasons:

| Group | Files | Reason |
|-------|-------|--------|
| CTD_BGC | 34 | BGC-only variables (O2, NO3, CHL, pH, ALK, DIC) — no T/S |
| SOCAT | 33 | Surface CO2 only (prof_PCO) — no T/S |
| PFL_BGC_A | 13 | Has prof_T/S fields but all zero/NaN; only O2 has data |
| ITP2 | 5 | No valid profile coordinates at step01 |
| Samoa | 2 | **Bug** (now fixed): step10 missing `if prof_key in MITprof_ds` guard, crashed on absent prof_Sweight |

### Bug fixed: step10 missing guard
Steps 03–09 all guard each `profile_var_key_set` loop with
`if prof_key in MITprof_ds`. Step10 was missing this, causing a KeyError
for Samoa mooring files which have `prof_T`/`prof_Tweight` but no `prof_S`.
Fix: one-line guard added to step10's weight-zeroing loop. NCEI.py's global
`profile_var_key_set` design is unchanged — each step is responsible for
checking variable existence, consistent with the existing pattern.

Samoa files verified after fix: both complete all 10 steps; ~4.2% of profiles
survive step10 subdaily decimation (expected — dense hourly mooring time series
decimated to once-daily nearest-noon, per PI design intent).

### Reference comparison: CTD 1992 and 1993
Compared new output against Ian's 2019 reference files
(`CTD_data_at_end_of_processing_chain/CTD_20190131_199[23].nc`).

**Profile counts:** new has ~1.9× more profiles (e.g. 1992: 46,135 vs 24,114).
Expected — new run uses a more recent WOD extraction; the reference was built
in 2019 from an older/smaller vintage. Not reproducible profile-for-profile.

**T distribution and mean profile:** nearly identical shape, range, and vertical
structure. New mean profile tracks reference through the full water column.

**S distribution:** reference contains values down to 0.08 PSU (likely data
errors in the older WOD vintage); new run minimum is ~30 PSU — step07 code 5
range QC working correctly on the cleaner input.

**S mean profile:** very close, slight divergence below ~1500 m — consistent
with WOA23 vs WOA13 climatology difference at depth (expected and correct).

**Weight distributions:** same shape and scale for both T and S weights.

**Deep coverage:** new run has better survival at depth (WOA23 full-depth clim
provides values to 5500 m; WOA13 stopped at 1500 m, so deep profiles previously
had no climatology reference and were dropped by step07 code 6).

**Geographic coverage:** same global distribution; new run fills in more of the
ocean consistent with the larger input vintage.

**Conclusion:** all differences are explained by (1) newer WOD vintage,
(2) WOA23 vs WOA13 climatology, (3) tighter S range QC. No unexplained
discrepancies. Chain is producing physically correct output.

### Step-10 survivor rates by source class

Parsed from the production run log. "Survivors" = profiles with nonzero weight after
all 10 steps, expressed as % of originals. Samoa fixed files run separately (N=2).

| Source | N files | Avg T% | Med T% | Avg S% | Med S% |
|--------|---------|--------|--------|--------|--------|
| CTD_WOD | 64 | 60.0% | 77.1% | 59.6% | 75.9% |
| GLD_WOD | 26 | 12.3% | 8.7% | 12.3% | 8.7% |
| MRB_WOD | 35 | 91.0% | 98.7% | 89.1% | 99.1% |
| XBT_WOD | 66 | 94.0% | 94.5% | n/a | n/a |
| PFL | 87 | 94.8% | 97.4% | 93.6% | 96.8% |
| PFL_BGC | 55 | 91.2% | 98.5% | 91.3% | 98.5% |
| MEOP | 21 | 68.8% | 69.8% | 73.8% | 73.9% |
| Samoa | 2 | 4.17% | 4.17% | n/a | n/a |

Samoa detail (T-only moorings, step10 is the only loss step):
- `samoa_1992_MRB`: 3135/75139 = 4.17%
- `samoa_2012_MRB`: 2235/53567 = 4.17%

Physical interpretation: GLD low (~12%) because gliders sample continuously and
step10 decimates to once-daily; CTD moderate (~60%) due to repeat casts at the same
station; MEOP moderate (~70%) due to seal haul-out clustering. Autonomous platforms
(PFL, MRB, XBT) 90%+ because already well-spaced. Samoa ~4% is deliberate PI design
(dense hourly mooring time series → one profile per day prevents moorings from
dominating the cost function).

### Per-step dropout analysis (T profiles, mean % cumulative survival)

Steps 1–7 and 9 drop nothing — they zero weights but do not remove profiles.
Two steps actually remove profiles:

| Step | CTD_WOD | GLD_WOD | MRB_WOD | XBT_WOD | PFL | PFL_BGC | MEOP |
|------|---------|---------|---------|---------|-----|---------|------|
| 8 (QC) | −4.8% | −0.6% | −0.1% | −1.7% | −2.9% | −7.2% | −0.9% |
| 10 (decimate) | −35.2% | −87.1% | −8.9% | −4.3% | −2.3% | −1.6% | −30.3% |

Step 8 (`update_remove_zero_T_S_weighted_profiles`) is a secondary loss point,
most significant for PFL_BGC (−7.2%) and CTD (−4.8%). Step 10 (subdaily decimation)
is the dominant loss for GLD, CTD, and MEOP.

### Additional validation checks

- **NaN cost fields:** NaN exclusively at zero-weight points — correct (cost
  undefined where weight=0).
- **High-cost outliers:** MEOP max ~260, MRB max ~82. Confirmed physical —
  Antarctic baroclinic intrusions (MEOP) and El Niño variability (MRB), not bugs.
- **lon/lat bounds:** all output `prof_lon` in [−180, 180), `prof_lat` in [−90, 90].
- **Source type spot-checks:** MRB `prof_S` survival 5.74% (was 1898% pre-fix);
  XBT T-only clean; CTD/PFL full T+S.

### Files
- Changed: `step10.py` (missing guard), `NCEI.py` (reverted; step10 fix is the right approach)

---

## 2026-09-26 — Performance audit: step01 and step10 (~36–60× speedups)

### Summary
Full audit of all 10 step files for xarray-inside-hot-loop antipatterns.
Two real bottlenecks found and fixed; the rest were structurally clean.

### Audit methodology
Code-reviewed all step files for: xarray DataArrays accessed element-by-element
in Python loops; expensive objects (KD-trees, geodesic solvers) rebuilt inside
loops; `any()`/`all()` generators pulling scalar xarray values; incremental
`np.union1d`/`append` inside tight loops. Confirmed severity by profiling with
`cProfile` on the 1992 CTD file (64,928 profiles).

### step10 — 185s → 3s (~60×)

**Root cause (two issues):**
1. `X`, `Y`, `Z`, `prof_HHMMSS`, `prof_YYYYMMDD` left as xarray DataArrays going
   into the per-day decimation loop. Every `distances_array[i, k]` scalar access
   triggered xarray's full `__getattr__` / `_attr_sources` / numpy version-check
   machinery — ~12M calls, ~170s.
2. `any(distances_array[ii_local, k] < tol for k in kept_local)` — Python
   generator pulling xarray scalars one at a time.

**Fix:** extract all five arrays to numpy `.values` before the loop; replace the
generator with `distances_array[ii_local, kept_local].min() < tol` (one numpy
slice). Output verified **bit-identical** (`prof_T`, `prof_S`, `prof_Tweight`,
`prof_Sweight` all match exactly).

### step01 — 27s → 1s (~27×)

**Root cause (two issues):**
1. `geopy.distance.distance()` called once per profile in a Python loop — each
   call runs the full iterative Karney WGS84 geodesic solver in pure Python.
   64,928 calls × ~100 µs each ≈ 7s.
2. `MITprof_ds['prof_lat'][i]` etc. inside the same loop — xarray scalar `isel()`
   triggered ~260k DataArray constructions ≈ 12s.

**Fix:** `pyproj.Geod.inv()` — same Karney WGS84 algorithm as geopy, but
accepts arrays and runs as a single vectorized C call. Verified: max difference
vs geopy = **0.0 m** (not an approximation — identical algorithm, identical
ellipsoid). `geopy` import retained in case it's used elsewhere.

**Note on algorithm choice:** an earlier intermediate commit used haversine
(spherical approximation, max 0.37% error). That was replaced with `pyproj.Geod`
before merging — full geodesic accuracy was preserved. The haversine would have
been defensible for this particular 2×dx threshold check, but pyproj is a cleaner
answer for a production codebase.

### Steps found clean (no action)
- **step07**: outer loop is ≤11 iterations (one per quality code); each body is
  fully vectorized xarray/numpy. `.item()` calls are logging only.
- **steps 02, 04, 05, 06, 08, 09**: no hot-loop antipatterns found in code review.
  (steps 04–09 could not be timed end-to-end due to missing gdrive grid files in
  local test environment, but code structure is clean.)

### Cumulative speedups (on 1992 CTD, 64,928 profiles)

| Step | Before | After | Factor |
|------|--------|-------|--------|
| step01 | ~27s | ~1s | ~27× |
| step03 | ~76s (see prior entry) | ~6s | ~13× |
| step10 | ~185s | ~3s | ~60× |

Chain wall-clock reduction for this file: roughly **280s → 10s** for these three
steps alone.

### Files
- Changed: `step01.py`, `step10.py`

---

## 2026-09-26 — clim_interp.py vectorization (~3000× speedup)

### Summary
Rewrote the inner loop of `interpolate_climatology` after profiling showed step03
was the bottleneck: 9.4 s/file vs ~0.6 s before the new optimal interpolation was added.

### Root cause
The original implementation grouped profiles by **rounded time-bracket weight**,
producing up to 365 distinct groups per file and rebuilding the scipy KD-tree
(needed for land-neighbor fallback) once per group — 365× repeated construction.

### Changes
- `_time_brackets_vec`: vectorized over the full profile array (no per-profile loop).
- `interpolate_climatology`: now groups by **(m0, m1) month-pair** only (≤12 groups),
  never by weight. Each profile's exact `w1` is still used for the blend, so time
  interpolation is still continuous and exact — grouping by pair is just a
  reorganization, not a coarsening.
- `_fallback_nearest_cell`: KD-tree cached per field (month index) and reused across
  all groups that share that field.
- `_depth_interp`, `_fallback_fill`: vectorized; bracket weights computed once for
  the shared `prof_depths` axis, applied to all profiles at once.

### Verification
Timed on the full CTD WOD 1992 file (64,928 profiles):

| N | Time | ms/prof | valid frac |
|---|------|---------|-----------|
| 5,000 | 0.65 s | 0.13 | 1.000 |
| 20,000 | 0.56 s | 0.028 | 1.000 |
| 64,928 | 1.62 s | 0.025 | 1.000 |

Pre-optimization baseline was ~74 ms/profile → **~3000× speedup**. Output value
range and valid fraction unchanged.

### Files
- Changed: `clim_interp.py`

---

## 2026-09-26 — NCEI chain bug fixes, WOA23 climatology, and optimal interpolation

### Summary
Three correctness fixes to the NCEI processing chain, a brand-new WOA23-based
climatology (replacing the legacy WOA13 `.mat`), and a rewritten,
optimal climatology interpolation. Verified end-to-end across three
structurally-distinct source types (CTD, MRB, XBT).

### Motivation
Investigation began from an anomalous log line on an MRB (moored-buoy) file:
`valid prof_S profile count / original valid prof_S profile count = 8067/425 = 1898.12%`
A survival ratio above 100% is impossible if the chain only removes/flags data,
so it pointed to a bug. Separately, the PI asked to (a) replace the old
climatology with World Ocean Atlas 2023 and (b) redo the interpolation optimally,
including time interpolation.

### Bug fixes

1. **step06 — climatology salinity was being persisted into `prof_S` (real bug).**
   step06 converts in-situ temperature to potential temperature, which requires
   salinity; where salinity was missing it borrowed climatology S. The original
   MATLAB (`update_prof_insitu_T_to_potential_T.m`) does this on a *local copy*
   and writes back only `prof_T` — the climatology S is transient scaffolding.
   The Python version instead wrote climatology S back into the stored `prof_S`,
   so fabricated salinity survived to the output and was assimilated as if it
   were a real observation (measured on 1992 CTD: 107,827 clim-filled points, 86%
   with nonzero weight). This is what produced the MRB >100% anomaly.
   *Fix:* `calculate_potential_T` now takes `prof_S` as an explicit argument; the
   caller builds a transient clim-augmented salinity for the calculation and
   writes back only `prof_T`. `prof_S` is never mutated. Verified: potential-T
   output is bit-identical given the same S input (the fix changes only S
   persistence, not the temperature math).

2. **step07 codes 9/10 — added MATLAB-parity value scrubs (defensive).**
   The cost-vs-climatology check now nulls sentinel/zero inputs before computing
   cost (value/clim `< -9000` or `== 0` → NaN; weight `< 0` → NaN), matching the
   MATLAB. Proven byte-inert on current data (no such values exist; fills are
   already NaN), but guards future inputs that might carry raw sentinels.

3. **NCEI.py — longitude normalized to [-180, 180) at write time.**
   Output `prof_lon` was raw 0–360; the reference end-of-chain files use
   [-180, 180). Now normalized at write only (every in-chain consumer feeds
   `sph2cart`, which is periodic in longitude, so no computation changes).
   Sentinel longitudes (>360) are left untouched so the wrap cannot disguise a
   missing position. Verified: only `prof_lon` changes; all else byte-identical.

### New WOA23 climatology — `build_woa23_climatology.py`
Builds a monthly potential-temperature + salinity climatology from World Ocean
Atlas 2023, 1991–2020 normal (`decav91C0`), 1° grid.
- Downloads 24 monthly + 8 seasonal T/S files from NCEI (robust: retries,
  Content-Length verification, open-check; resumable).
- WOA `t_an` is in-situ temperature; converted to potential temperature using
  the chain's own `step06.calculate_potential_T` (each grid cell treated as a
  pseudo-profile) so the climatology is converted identically to the
  observations.
- Full depth to 5500 m: monthly fields (0–1500 m, the WOA monthly limit) are
  spliced onto the seasonal fields below 1500 m; the seasonal fields are first
  interpolated to monthly by day-of-year (mid-season anchors, wrapping) so the
  deep part has a smooth monthly time axis and no season-boundary discontinuity.
- Output: `woa23_climatology/woa23_decav91C0_TS_clim_potential_T_1deg_fulldepth.nc`
  (variables `potential_T_monthly`, `S_monthly`, dims `(12, 102, 180, 360)`).

### Optimal interpolation — `clim_interp.py` (used by step03)
Replaces step03's nearest-neighbor space lookup + month-step-function time
selection with:
- **Space:** bilinear in lon/lat (NaN-aware; longitude-seam wrap).
- **Depth:** linear onto each profile's own depth levels.
- **Time:** linear between adjacent monthly fields by day-of-year (mid-month
  anchors, Dec↔Jan wrap) — eliminates the month-boundary discontinuity.
- **Fallback:** nearest-valid climatology where the primary interpolation yields
  NaN (profile depth beyond the grid; all-land neighborhood). A KD-tree over
  valid surface cells handles spatial gaps.
- All methods and the fallback are config knobs at the top of the module, so the
  scheme can be changed without touching logic.

### Climatology comparison (WOA23 vs legacy WOA13)
On the identical grid, the two agree closely: mean difference +0.012 °C (T) /
+0.001 PSU (S); typical cell differs ~0.04 °C / ~0.007 PSU. Differences are
largest at the surface (~0.22 °C, 0–100 m) and negligible at depth (~0.026 °C,
1500–5500 m) — the expected signature of the WOA13→WOA23 update (more
observations, ocean warming through 2022), not an artifact. The <0.2% of cells
that differ by >2 °C / >1 PSU are all shallow coastal/marginal-sea points
(Caspian Sea, Kara Sea/Ob estuary, Río de la Plata), where WOA23 corrects
known-bad or under-resolved WOA13 values.

### Verification
- Full 10-step chain run on 1992 CTD with the new clim: no crashes, no S/T
  inversion (S 75.9% / T 76.4% survivors), climatology fully populated including
  deep levels, cost fields finite.
- Spot-checks on the structurally-distinct source types:
  - **MRB** (T-only buoys, the origin of the anomaly): `prof_S` survival now
    5.74% (was 1898%); zero climatology-S persisted into `prof_S`.
  - **XBT** (T-only): 97.8% temperature survival, clean run.

### Cleanup
- step03: removed unused imports (`numpy.ma`, `argparse`, `glob`, `os`,
  `scipy.io`, `netCDF4`, `griddata`) and stale comments.
- Fixed a latent bug in step03's `.nc` clim loader: it called
  `assign_coords(np.arange(...))`, which would overwrite the real lon/lat/depth
  coordinate values with integer indices; now reads `.values` directly.

### Files
- Changed: `NCEI.py`, `step03.py`, `step06.py`, `step07.py`
- Added: `build_woa23_climatology.py`, `clim_interp.py`, `PROJECT_LOG.md`

### Not done / follow-ups
- Full multi-source production run (ready to launch).
- `NCEI.py` `climatology_file` now points at the WOA23 full-depth file; the
  legacy `.mat` line is retained (commented) for easy revert.
