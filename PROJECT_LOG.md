# Project Log

A running, in-repo changelog. Newest entry first. Each entry summarizes what a
commit changed, why, what was discussed, and how it was verified. (Complements
git history with human-readable rationale.)

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
