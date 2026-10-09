# Handoff Context — ECCO-Insitu-Python

**Branch:** `devel_v3`  
**Last commit:** `3c12a6b`  
**Date:** 2026-10-09

---

## Current state

### NCEI processing chain

- The NCEI 10-step chain runs via `run_ncei.sh`, followed by post-chain steps: **step 11 spatial scaling** (on by default, `-S` to skip), diagnostic plots (`-p`), **step 12 compaction** (off by default, `-c`).

### Step 11 spatial scaling (changed 2026-10-09)

- `step11_spatial_scaling.py` now scales **T and S independently**: n_V = profiles with any nonzero `prof_Vweight` per (10242-bin, month) across all sources; `prof_Vweight *= 1/n_V`. Writes `prof_Tspatial_scaling_factor` / `prof_Sspatial_scaling_factor`. Handles T-only files (Samoa). Refuses to run on already-scaled files.
- This is a best guess at what the team wants; may change with new instructions. Variables scaled are set by `SCALED_VARS` (U/V in Samoa are not scaled).
- Old combined-count version (matches original MATLAB) kept as `step11_spatial_scaling_combined_LEGACY.py`, not in the chain.
- `profile_files_NCEI_processed_20261007_144529` was scaled with the **LEGACY** version, and its two Samoa files were skipped (no `prof_Sweight`) — they are still unscaled step-10 files. Re-run new step 11 on a copy of `..._20261007_144529_step10_backup`.
- `plot_spatial_scaling.py` is per-variable (28 PNGs, `_T`/`_S` tags), reusing step 11's counting.
- Step 11 only changes weights; T/S/clim/bins are untouched, so the globe/scatter apps and geodesic binning plots are unaffected. `prof_Tcost`/`prof_Scost` are from step 07 (pre-scaling weights) and are not recomputed.
- Most recent completed run: `profile_files_NCEI_processed_20261006_192019` (in `/Users/brucel/ecco/yip/processed_by_NCEI_profile_data/`)
- **Known issue with this run**: the `booL_mask_notnull_name` typo in `tools.py` (capital L) caused step 08 crashes on some GLD files. This was fixed in commit `3c12a6b` but the run was not restarted after the fix — so some GLD files from this run may be missing output.
- **`iINTERP` dimension**: Fortran code needs `prof_interp_*` variables to have dims `(iPROF, iINTERP)` with `iINTERP=1`. Fixed in `tools.py:MITprof_write_to_nc` — all future runs will produce the correct shape.

### Timeseries / visualization data

Two timeseries stores exist under `/Users/brucel/ecco/yip/profile_data/z_profile_file_analysis/`:

1. **`allsource_bin_timeseries_20260930_164848/`** — new layout (directory contains zarr + json + will need profiles.db). Built from the `20260930_164848` NCEI run. Compatible with `app_allsource_globe.py`. **Missing `_profiles.db`** — SQLite profile records were added to the build scripts after this was built, so bin-click source tables won't work until rebuilt.

2. **`allsource_bin_timeseries_with_anomalies_20260930_164848_metadata.json` + `.zarr`** — old flat layout (files loose in `z_profile_file_analysis/`, not inside a container directory). Built before the directory-structure change. **Not compatible with `app_allsource_anomaly_globe.py`** as-is — the app now expects files inside a `<stem>_YYYYMMDD_HHMMSS/` directory.

### What needs to happen next

1. **Rebuild both timeseries** from the `20261006_192019` NCEI run (or a cleaner one) to get:
   - new directory layout
   - SQLite `_profiles.db` included
   - anomaly arrays in the same store
   - Data from the fixed chain (iINTERP + booL typo fixed)

   Commands (from repo root):
   ```bash
   python profile_analysis/build_allsource_bin_timeseries.py \
       --source_root /Users/brucel/ecco/yip/processed_by_NCEI_profile_data/profile_files_NCEI_processed_20261006_192019

   python profile_analysis/build_allsource_bin_timeseries_with_anomalies.py \
       --source_root /Users/brucel/ecco/yip/processed_by_NCEI_profile_data/profile_files_NCEI_processed_20261006_192019
   ```

2. **Test `app_allsource_anomaly_globe.py`** — this has not been tested end-to-end with the new data layout yet.

---

## Pending / shelved issues

### MRB_WOD data ambiguity (awaiting team feedback)
Under `/Users/brucel/ecco/yip/profile_data/Interp_Profiles/MRB_WOD/` there are three overlapping datasets:
- `MRB_WOD/` (top-level): 1.07M profiles across all years
- `MRB_WOD/36_depths/`: 44K profiles/year on 36-level grid
- `MRB_WOD/97_depths/`: same 44K profiles as 36_depths, interpolated to 97-level grid

The top-level and `97_depths` sets overlap ~86% by location/date, but T values differ by up to ~0.8°C for matched profiles — they are different versions of the data, not just different grids. The chain currently ingests all three via `rglob`, triple-counting many profiles. Decision needed on which folder(s) to use before the next production run.

### Globe rotation gimbal lock (deferred)
When panning near the poles (occasionally near other continents), the globe zoom briefly snaps to a higher level then recovers. Root cause: gimbal lock in Euler-angle rotation math of the clientside drag handler. Fix: rewrite clientside rotation using quaternions. Deferred.

### Flat map zoom-out jerkiness (deferred)
When zooming out on the flat map, plot limits jump rather than updating smoothly. Root cause: Plotly redraws the entire choropleth figure on bounds change — no incremental update path. Possible mitigation: snap to discrete zoom levels. Deferred.

---

## Key file locations

| What | Path |
|------|------|
| NCEI chain controller | `NCEI.py` |
| Run script | `run_ncei.sh` |
| Decimation (step10) | `step10.py` |
| Spatial scaling (step11) | `step11_spatial_scaling.py` (+ `_combined_LEGACY.py`) |
| Scaling diagnostic plots | `plot_spatial_scaling.py` |
| Compaction (step12) | `step12_profiles_compact.py` |
| Build timeseries | `profile_analysis/build_allsource_bin_timeseries.py` |
| Build timeseries w/ anomalies | `profile_analysis/build_allsource_bin_timeseries_with_anomalies.py` |
| Globe app | `profile_analysis/app_allsource_globe.py` |
| Anomaly globe app | `profile_analysis/app_allsource_anomaly_globe.py` |
| WOA23 climatology | `/Users/brucel/ecco/yip/woa23_climatology/woa23_decav91C0_TS_clim_potential_T_1deg_fulldepth.nc` |
| NCEI processed output | `/Users/brucel/ecco/yip/processed_by_NCEI_profile_data/` |
| Timeseries stores | `/Users/brucel/ecco/yip/profile_data/z_profile_file_analysis/` |
| Input profiles | `/Users/brucel/ecco/yip/profile_data/Interp_Profiles/` |

---

## Architecture notes

- All sources share a 97-level reference depth grid; shallower instruments are NaN-padded.
- Zarr layout: 2D `(n_months, 97)` arrays per source per bin, with shared `ym_index` (YYYYMM int32) for alignment.
- Profile records (file, prof_idx, date, lon, lat, qual) are stored in SQLite (`_profiles.db`), not in the JSON, to keep startup memory manageable.
- Apps auto-select the most recent timestamped store directory; `--data <path>` overrides.
- `aggregate_bins` uses `ProcessPoolExecutor(fork)` — known risk: fork+OpenMP deadlock. Fix: `OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 python build_allsource...`
