# Project Log

A running, in-repo changelog. Newest entry first. Each entry summarizes what a
commit changed, why, what was discussed, and how it was verified. (Complements
git history with human-readable rationale.)

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
