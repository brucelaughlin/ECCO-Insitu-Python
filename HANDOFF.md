# Handoff Context — ECCO-Insitu-Python

**Branch:** `devel_v3`
**Last commit:** `1aa5bd1` (much of the work below is not yet committed — see the end)
**Date:** 2026-10-12

---

## Current state

### Inputs (`/Users/brucel/ecco/yip/profile_data/Interp_Profiles/`)
246 files, all on one 97-level depth grid (2–6000 m). Deliberately moved to
`profile_data/Raw_Profile_Data_Unused_By_Ecco/`: Argo A and R files (Argo is
delayed-mode **D only**), PFL_BGC (95% duplicate of PFL T/S), CTD_WOD `*_BGC` (no
T/S), SOCAT (no T/S), MRB_WOD `36_depths/` and `97_depths/` (overlapping
duplicates), Samoa (own deep grid). Two ITP L2 files (2008, 2018) are empty and
three (2012, 2014, 2017) have bad salinity; the chain rejects all five.

### Chain (`run_ncei.sh` → `NCEI.py`, then post-chain steps)
| Step | What | Notes |
|---|---|---|
| 01–04 | grid/tile points, geodesic bins, WOA23 climatology, σ | step 01 assigns LLC90 cells and 30×30 tiles |
| 05 | potential T | |
| 06 | QC (zero weights) | runs **before** γ since 2026-10-12, so its cost tests use (ΔV/σ)²; skips the cost tests poleward of 60° |
| 07 | drop all-zero-weight profiles | |
| 08 | trim empty depth levels | **off** by default (`remove_extraneous_depth_levels = False`) → every file on the 97-level grid |
| 09 | decimate sub-daily repeats | removes ~92% of glider and ~61% of ITP profiles |
| 10 | γ (grid-cell area factor) | moved here from step 05 on 2026-10-12 |
| 11 | spatial scaling, T and S separately | rewrites files in parallel (8 workers, `-n`; ~3× faster); `-S` to skip |
| 12 | compaction | off by default (`-c`) |
| 13 | ECCO-ready profile files | `-A` to skip; see below |

Output: `processed_by_NCEI_profile_data/profile_files_NCEI_processed_<ts>/` with
`yearly/` (per source and year), `ecco_arrays/` (step 13 model-ready files only)
and `documentation/` (step 13 manifest, trade-off sweep, `PROFILES_SIZE.h` values
and `data.profiles` list).

High-latitude uncertainties (2026-10-14), set in `NCEI.py`, step 04:
- **North of 60°N:** σ from `profile_data/high_latitude_uncertanties/` (`Theta_rmse_se`,
  `S_rmse_se`; rmse confirmed by the team). EASE2 100 km, Arctic only (59.4–89.4°N);
  nearest cell within 300 km.
- **South of 60°S:** no file yet; for now (team instruction) a constant σ_T = 10 °C,
  σ_S = 5 psu at all depths (`southern_sigma_const_dict`).
- `prof_Vuncertainty_source`: 0 = CTD σ, 1 = Arctic file, 2 = southern constant. Step 06
  runs its cost tests on 1 and 2 (`clim_cost_qc_with_high_lat_sigma`); with the southern
  σ that large, they reject essentially nothing there (MEOP 2010: none).
- Effect: Arctic S σ often larger at high-misfit points (fewer large costs; the QC removes
  the rest). Southern weights drop to ~1/100 (T) and ~1/25 (S) of a σ = 1 profile.
- Not yet in a full run with the southern constant (the 2026-10-14 10:44 run has the Arctic σ only).

Date/time clean-up (in `NCEI.py`, before step 01): out-of-range `HHMMSS` /
`YYYYMMDD` values from WOD rounding (seconds or minutes = 60, 24:00:00, day past
month end) are rolled forward with real date arithmetic, so they agree with
`prof_date` (≈3,900 profiles, mostly gliders).

### Step 13 / `optimize_profile_files.py` — model-ready files
Target: ECCO v4r6 (LLC90, 30×30 tiles; 117 tiles, 113 active). The model
reserves per tile `NFILESPROFMAX` × `NOBSGLOB` profile slots (× levels), where
`NOBSGLOB` is the most profiles any one file puts in any one tile. The program
minimizes that product: one provenance per file (ITP L2/L3 separate, XBT flag +
noflag together), chronological, files labelled in whole calendar years (several
whole years, or parts of one too-dense year with the year's label). Default
`NOBSGLOB = 'auto'` = fewest files within 5% of the minimum product. Needs only
the LLC size and tile size — not blank tiles or bathymetry (the model ignores
profiles on land / in blank tiles itself). Writes the `PROFILES_SIZE.h` values and
the `data.profiles` file list. **Plain-language guide:
`docs/profile_file_optimization.md`.** Validate with
`diagnostics/check_ecco_arrays.py <run>`.

Result for the 2026-10-12 data: 18 files, `NOBSGLOB` 43,469 (product 782,442;
minimum 762,750 at 27 files; lower bound 599,691). Argo gets 3 files
(1997–2013, 2014–2021, 2022–2026), close to the colleague's own split.

### Runs
- **`profile_files_NCEI_processed_20261012_172848`** — newest; new step order,
  new layout. `ecco_arrays/` re-made with the current step 13 on 2026-10-12
  (summary files moved to `documentation/`):
  **18 files, set `NOBSGLOB = 43351`, `NFILESPROFMAX = 18`** (5,941,130
  profiles; `check_ecco_arrays.py` PASS, log `z_logs/check_ecco_arrays_20261012_172848.log`).
  `check_run_files.py` not yet run on it.
- `..._20261012_103552` — old step order (γ before QC), old flat layout; passed
  `check_run_files.py`. Used for the cost reports.
- `..._20261009_151941` — first run from the cleaned-up inputs; old step order.

### Diagnostics (`diagnostics/`, see `diagnostics/README.md`)
`check_run_files.py` (structure, ~15 s), `check_run_qc.py` (science QC, ~1 h —
run with `nohup`), `check_ecco_arrays.py`, `report_profile_counts.py` (UCSD vs
processed counts, maps), `report_costs.py` (cost vs climatology with/without γ and
step 11 factors, sf_T vs sf_S, decimation), `report_spatial_scaling.py`.
All accept either run layout and read each run's step order from its log.

### Findings worth remembering
- **Large costs**: the extreme ones (up to ~10⁵) are all poleward of 60°, where
  step 06 skips its cost tests — bad climatology (Laptev Sea S = 25.7 at all
  depths) and bad data (MEOP 2019 Ross Sea T ≈ 30 °C, a 2003 Ross Sea XBT). With
  γ before QC (old order) ~27,600 mid-latitude profiles also passed QC with
  (ΔV/σ)² > 100; the new order removes that effect.
- **sf_T ≠ sf_S** for ~30% of profiles, because XBTs, moorings and seal tags add
  T-only profiles to bin-months.
- **446 profiles per run have no geodesic bin** (nearest LLC90 cell is inland or
  ice-shelf land); kept for now — see the to-do below.
- **4,841 profiles (0.08%)** sit exactly equidistant between two LLC90 cells at
  facet seams; step 01's choice is self-consistent, so no action.

---

## Pending

1. **Check the yearly files** of `..._20261012_172848`:
   `python diagnostics/check_run_files.py <run> --input_dir …/Interp_Profiles --log z_logs/ncei_log_20261012_172848.log`
   (its `ecco_arrays/` already passed `check_ecco_arrays.py`).
2. **Confirm with the model team**: the `PROFILES_SIZE.h` values step 13 prints;
   `NLEVELMAX` could drop from 140 to 97 (our depth grid), shrinking every tile's
   arrays by ~30%.
3. **To-do next time — drop profiles with no geodesic bin** (approved idea, not
   applied). Tested patch: `patches/step06_code12_drop_unbinned_profiles.patch`
   adds QC code 12, toggled by `zero_weights_of_unbinned_profiles` in `NCEI.py`.
   `git apply --check patches/step06_code12_drop_unbinned_profiles.patch`, then
   `git apply …`.
4. **Possible QC change for the team**: switch the cost tests back on poleward of
   60° with looser thresholds (would catch the bad seal / XBT data); the bad
   Laptev climatology needs its own fix.
5. **Timeseries stores for the globe apps** need rebuilding from a current run
   (`profile_analysis/build_allsource_bin_timeseries*.py --source_root <run>`;
   `yearly/` is found automatically).
6. Deferred app issues: globe rotation gimbal lock near the poles (fix: quaternion
   rotation); flat-map zoom-out jerkiness.

---

## Uncommitted work (as of 2026-10-12)
`NCEI.py`, `step01`–`step10` (renumbering, step 08 switch), `run_ncei.sh` (layout,
step 13), `optimize_profile_files.py`, `step13_pack_ecco_arrays.py`,
`diagnostics/` (all), `docs/`, `patches/`, `README.md`, this file, and earlier
edits to `csv_to_nc.py`, `tools.py`, the `profile_analysis/` apps and build
scripts.

---

## Key file locations

| What | Path |
|------|------|
| Chain controller / run script | `NCEI.py`, `run_ncei.sh` |
| Steps | `step01.py` … `step10.py`, `step11_spatial_scaling.py`, `step12_profiles_compact.py`, `step13_pack_ecco_arrays.py` |
| File-packing program | `optimize_profile_files.py` (guide: `docs/profile_file_optimization.md`) |
| Checks & reports | `diagnostics/` (see `diagnostics/README.md`) |
| Saved, tested patches | `patches/` |
| Globe apps and their build scripts | `profile_analysis/` |
| Input profiles | `/Users/brucel/ecco/yip/profile_data/Interp_Profiles/` |
| Unused inputs | `/Users/brucel/ecco/yip/profile_data/Raw_Profile_Data_Unused_By_Ecco/` |
| Chain output | `/Users/brucel/ecco/yip/processed_by_NCEI_profile_data/` |
| Reports | `/Users/brucel/ecco/yip/cost_reports/`, `/Users/brucel/ecco/yip/insitu_reports/` |
| Logs | `z_logs/` |
| WOA23 climatology | `/Users/brucel/ecco/yip/woa23_climatology/` |
| ECCO v4r6 config (tiles, PROFILES_SIZE.h, data.profiles) | github.com/ECCO-GROUP/ECCO-v4-Configurations, `ECCOv4 Release 6/` |

## Architecture notes (globe apps)
- Zarr layout: 2-D `(n_months, 97)` arrays per source per bin, with a shared
  `ym_index` (YYYYMM int32).
- Profile records are in SQLite (`_profiles.db`), not the JSON, to keep start-up
  memory manageable.
- Apps auto-select the most recent timestamped store; `--data <path>` overrides.
- `aggregate_bins` uses `ProcessPoolExecutor(fork)`; if it deadlocks, run with
  `OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1`.
