# ECCO In-Situ Processing Pipeline

Authors: Sweet Zhang, Ian Fenty, Bruce Laughlin  
Originally: Sweet Zhang, 2023-05-21  
Updated: 2026-10

---

## Overview

This repository contains a Python pipeline for post-processing in-situ ocean
profile data (from WOD, Argo, MEOP, mooring, and other sources) into a form
suitable for assimilation into the ECCO ocean state estimate. The pipeline:

1. Converts raw WOD CSV exports to NetCDF (`csv_to_nc.py`)
2. Runs a 10-step quality-control and preprocessing chain on each profile file (`NCEI.py`)
3. Applies spatial density scaling to the weights (step 11)
4. Packs the result into model-ready profile files that minimize the model's
   memory (step 13, `optimize_profile_files.py`; plain-language explanation in
   [`docs/profile_file_optimization.md`](docs/profile_file_optimization.md))
5. Checks and reports on the output (`diagnostics/`, see `diagnostics/README.md`)
6. Aggregates the processed profiles into geodesic bins for visualization, with an
   interactive globe app (`profile_analysis/`)

A support data package (LLC90 grid, WOA23 climatology, sigma files, geodesic
bin files) is required and is available on the project Google Drive:
https://drive.google.com/drive/folders/17h0qMS7vVimet8FXieGP1mWhnqnY0ljr

---

## Repository structure

```
ECCO-Insitu-Python/
  NCEI.py                         — main pipeline controller
  run_ncei.sh                     — shell wrapper for NCEI.py (recommended entry point)
  csv_to_nc.py                    — convert WOD CSV exports to NetCDF
  prebake_woa23_climatology.py    — offline climatology pre-interpolation (run once)
  clim_interp.py                  — climatology interpolation library
  high_latitude_sigma.py          — high-latitude T/S uncertainties for step 04
  step01.py … step10.py           — individual pipeline steps
  step11_spatial_scaling.py       — post-chain spatial density scaling (per T/S)
  step12_profiles_compact.py      — optional post-chain compaction
  step13_pack_ecco_arrays.py      — packs yearly output into ECCO-ready profile files (chain settings)
  optimize_profile_files.py       — the packing program: any LLC grid / tile size; trade-off report
  diagnostics/                    — checks and reports on chain output (validation, QC, profile counts,
                                    costs, step 11 figures); see diagnostics/README.md
  tools.py                        — shared utilities
  docs/profile_file_optimization.md — plain-language guide to the step 13 file packing (also .pdf)
  patches/                        — saved, tested changes not yet applied (see HANDOFF.md)
  profile_analysis/
    app_profile_map.py            — interactive map app (World / Arctic / Southern Ocean; bins or profiles)
    map_common.py                 — its projections, bin cells and zoom handling
    build_allsource_bin_timeseries_with_anomalies.py — builds the data store the app reads
    README_maps.md                — map app documentation
```

---

## Step 0 — Preprocessing: `csv_to_nc.py`

Converts a WOD `.csv` export to a set of per-year NetCDF files, with basic
validation of flags, coordinates, and units. Profiles failing validation are
excluded and logged.

```bash
python csv_to_nc.py -i <input_csv_dir> -d <output_dir>
```

Output filenames: `[stem]_[year].nc`

---

## Step 1 — NCEI processing chain: `NCEI.py` / `run_ncei.sh`

Processes each input NetCDF file through 10 sequential steps:

| Step | Function | Description |
|------|----------|-------------|
| 01 | `update_prof_and_tile_points` | Interpolates profiles onto LLC90 grid; assigns tile/grid indices |
| 02 | `update_spatial_bin_index` | Assigns geodesic bin IDs (`prof_bin_id_a` at 10242-bin, `prof_bin_id_b` at 2562-bin) |
| 03 | `update_monthly_mean_TS_clim` | Interpolates WOA23 T/S climatology to each profile's location, depth, and day-of-year |
| 04 | `update_sigmaTS` | Assigns T/S uncertainty σ (weight = 1/σ²) from the CTD sigma files; north of 60°N from the high-latitude uncertainty files (`high_latitude_sigma.py`), south of 60°S a constant σ_T = 10 °C, σ_S = 5 psu; records `prof_Vuncertainty_source` (0 CTD, 1 Arctic file, 2 southern constant) |
| 05 | `update_prof_insitu_T_to_potential_T` | Converts in-situ temperature to potential temperature; optionally fills missing S with climatology |
| 06 | `update_zero_weight_points` | Zeros weights on profiles failing quality criteria; optionally excludes high-latitude profiles from climatology cost. Runs before γ, so its cost tests use (ΔV/σ)² |
| 07 | `update_remove_zero_weighted_profiles` | Removes profiles with all-zero T and S weights |
| 08 | `update_remove_extraneous_depth_levels` | Removes trailing depth levels with no valid data (off by default: `remove_extraneous_depth_levels = False`) |
| 09 | `update_decimate_subdaily_profiles` | Decimates sub-daily sampling at the same location to once-daily |
| 10 | `update_gamma_factor` | Multiplies T/S weights by the grid-cell area factor γ (`prof_area_gamma`); moved after QC (was step 05) on 2026-10-12 |

### Post-chain steps (run by `run_ncei.sh` across the whole run directory)

| Step | Script | Default | Description |
|------|--------|---------|-------------|
| 11 | `step11_spatial_scaling.py` | on (`-S` to skip) | Spatial density scaling of T and S weights |
| — | `diagnostics/report_spatial_scaling.py` | off (`-p` to run) | Step 11 diagnostic figures |
| 12 | `step12_profiles_compact.py` | off (`-c` to run) | Rewrites files compactly in-place |
| 13 | `step13_pack_ecco_arrays.py` | on (`-A` to skip) | Packs the yearly files into ECCO-ready profile files (below) |

**Output layout** (inside each timestamped run folder):

```
profile_files_NCEI_processed_<YYYYMMDD_HHMMSS>/
  yearly/          one file per source and year (NCEI.py, steps 11, 12)
  ecco_arrays/     step 13: the model-ready profile files only
  documentation/   step 13 summary: ecco_arrays_manifest.json / .csv, optimize_sweep.csv,
                   PROFILES_SIZE_and_data.profiles.txt (NOBSGLOB / NFILESPROFMAX, data.profiles list)
```

**Step 11 — spatial density scaling.** Downweights densely sampled regions so
each (geodesic bin, window) contributes roughly one profile's worth of forcing.
For each variable V in `SCALED_VARS` (`T`, `S`), step 11 counts the profiles
with at least one nonzero `prof_Vweight` in each 10242-bin cell per window
(calendar month by default; `-w 10day` for 10-day windows), across all sources
combined, and multiplies each profile's `prof_Vweight` by 1/n_V. T and S are
counted separately, so T-only profiles (e.g. Samoa moorings) do not reduce the
S weights of other profiles in the same bin. Files carrying only some of the
weights are scaled for what they have. Adds `prof_Tspatial_scaling_factor` /
`prof_Sspatial_scaling_factor` and renames files to `__ncei_step_11.nc`.
Refuses to run on files that are already scaled — run it on a step-10 copy.
Rewriting the files is nearly all of its run time, so that part runs in parallel
(8 workers by default; `-n` to change).

The original MATLAB (`update_weights_based_on_spatial_density_TRYTOIMPLEMENT.m`)
counts all profiles together and applies one factor to both T and S (a Python
version of that, `step11_spatial_scaling_combined_LEGACY.py`, is in the git history
up to commit 1aa5bd1). Unlike the MATLAB, step 11 does not zero weights outside a
fixed range of scaling years.

```bash
python step11_spatial_scaling.py <run_dir> --dry_run   # preview
python step11_spatial_scaling.py <run_dir>             # apply in-place
python diagnostics/report_spatial_scaling.py <run_dir>  # 28 figures (14 per T/S)
```

Figures go to `<processed_root>/figures/figures_<YYYYMMDD_HHMMSS>/`, outside
the run directory.

**Step 13 — ECCO-ready profile files** (`optimize_profile_files.py`, called by
`step13_pack_ecco_arrays.py` with the chain's settings). New to this? Start with
[`docs/profile_file_optimization.md`](docs/profile_file_optimization.md). The MITgcm profiles package
allocates per-tile arrays of `NFILESPROFMAX` × `NOBSGLOB` (× levels), where
`NOBSGLOB` is the most profiles any one file puts in any one tile; the program
arranges the files to keep that product (i.e. the empty space) small: one
provenance per file (ITP split into L2 / L3), chronological, each file labeled in
whole calendar years (`YYYY0101`–`YYYY1231`; several whole years, or — if one year
is too dense — one part of a single year, e.g. `XBT_WOD_001_19920101_19921231`,
`XBT_WOD_002_19920101_19921231`). It sweeps `NOBSGLOB`, prints the trade-off
(files, product, % empty) and by default takes the fewest files within 5% of the
minimum (2026-10-12 data: 18 files, `NOBSGLOB` 43,469). It writes the
`NOBSGLOB` / `NFILESPROFMAX` to set and the `data.profiles` file list. Works for
any LLC grid size N and per-rank tile size (only N, the tile size and the grid's
cell coordinates are needed; which tiles a model run blanks doesn't matter):

```bash
python optimize_profile_files.py <profiles_dir>                         # report only
python optimize_profile_files.py <profiles_dir> --write <out_dir>       # write the files
python optimize_profile_files.py <dir> --nobsglob 35000 --year_rule none
python optimize_profile_files.py <dir> --llc 270 --tile 30 --grid_dir <llc270 grid dir>
```
Check with `python diagnostics/check_ecco_arrays.py <run_dir>`.

### Checking a run

Two quick checks should pass on every new run before its output is used (details
in [`diagnostics/README.md`](diagnostics/README.md)); both print PASS / FAIL and exit 0 on PASS:

- **`diagnostics/check_run_files.py`: are the yearly files complete and well-formed?**
  Every input has an output or a logged failure; all files at the same step and on
  one depth grid; weights finite, ≥ 0 and never on fill values; every profile used;
  bins and dates valid; and the step 11 factors sum to exactly 1 in every bin-month,
  for T and S. About 15 s.
- **`diagnostics/check_ecco_arrays.py`: are the model-ready files a faithful repacking?**
  Every yearly profile appears exactly once with identical data; one provenance
  per file, in time order, with whole-year labels; no tile over `NOBSGLOB`; the
  `PROFILES_SIZE.h` values in `documentation/` match the files. About a minute per GB.

```bash
RUN=/Users/brucel/ecco/yip/processed_by_NCEI_profile_data/profile_files_NCEI_processed_<ts>
python diagnostics/check_run_files.py $RUN --input_dir /Users/brucel/ecco/yip/profile_data/Interp_Profiles --log z_logs/ncei_log_<ts>.log
python diagnostics/check_ecco_arrays.py $RUN
```

A slower scientific check (`check_run_qc.py`, ~1 h) and the reports (profile
counts, costs, step 11 figures) are described in `diagnostics/README.md`.

### Required support files

All paths are set in `NCEI.py` under the `NEED PATHS / PARAMETERS` block:

| Variable | Description | Location in support package |
|----------|-------------|---------------------------|
| `grid_dir` | LLC90 grid directory | `grid_llc90/` |
| `sphere_bin_dir` | Geodesic bin assignment files | `grid_llc90/sphere_point_distribution/` |
| `climatology_file` | WOA23 full-depth T/S climatology (NetCDF) | `woa23_climatology/woa23_decav91C0_TS_clim_potential_T_1deg_fulldepth.nc` |
| `sigma_file_dict` | CTD T/S uncertainty fields (binary) | `CTD_sigma_TS/` |

Pre-baked climatology files (optional but strongly recommended for performance)
are picked up automatically: every `*_prebaked.nc` in `prebaked_clim_dir`
(set in `NCEI.py`) is loaded and matched to input files by the depth grid
stored in the file (`obs_depth`), so adding a new grid needs no code change.
Generate them once with:

```bash
python prebake_woa23_climatology.py \
    --input_dir <path_to_Interp_Profiles> \
    --source    <path_to_woa23_fulldepth.nc> \
    --out_dir   <path_to_woa23_climatology_dir>
```

This pre-interpolates the WOA23 climatology onto each distinct observation depth
grid found in the input files, eliminating per-profile vertical interpolation at
runtime (step 03). Without pre-baked files the chain still runs correctly but
is significantly slower.

### Running the chain

The recommended entry point is `run_ncei.sh`, which timestamps the output
directory and tees all output to a log file:

```bash
./run_ncei.sh                              # default input/output paths
./run_ncei.sh -i /path/to/input           # override input directory
./run_ncei.sh -n 8                        # use 8 parallel workers
./run_ncei.sh -i /path/to/input -n 8 -d /path/to/output
./run_ncei.sh -S                          # skip step 11 spatial scaling
./run_ncei.sh -w 10day                    # step 11 with 10-day windows
./run_ncei.sh -p                          # also make step 11 diagnostic figures
./run_ncei.sh -c                          # also run step 12 compaction
./run_ncei.sh -A                          # skip step 13 ECCO-ready packing
```

Parallelism default: `min(8, cpu_count - 1)` workers. Sequential mode: `-n 1`.
On a modern workstation, ~420 files complete in approximately 10–15 minutes
with 8 workers.

Alternatively, call `NCEI.py` directly:

```bash
python NCEI.py -i <input_dir> -d <output_dir> [-n <n_workers>]
```

### Output

Each run writes a timestamped folder containing `yearly/`, `ecco_arrays/` and
`documentation/` (see "Output layout" above). `yearly/` has one NetCDF file per input file, under
a sub-folder named by source (e.g. `yearly/CTD_WOD/`, `yearly/PFL/`), all on the
same 97-level depth grid. Each file contains the original profile data plus added
fields including:

- `prof_bin_id_a`, `prof_bin_id_b` — geodesic bin assignments
- `prof_Tclim`, `prof_Sclim` — WOA23 climatology at each profile
- `prof_Tweight`, `prof_Sweight` — quality weights (spatially scaled after step 11)
- `prof_Tspatial_scaling_factor`, `prof_Sspatial_scaling_factor` — step 11 factors
- `prof_Tcost`, `prof_Scost` — climatology misfit cost (computed in step 06, before γ (step 10) and step 11 are applied; not updated by them)

### Parameter reference

Parameters set in `NCEI.py` under the `NEED PATHS / PARAMETERS` block:

| Parameter | Default | Description |
|-----------|---------|-------------|
| `llcN` | 90 | LLC grid resolution (90 or 270) for step 01's grid / tile assignment |
| `wet_or_all` | 1 | Grid interpolation: 0 = nearest wet point, 1 = all points |
| `respect_existing_zero_weights` | False | Preserve pre-existing zero weights |
| `new_floor_dict` | T=0, S=0.005 | Minimum weight floors |
| `apply_gamma_factor` | True | Apply area-based gamma correction to sigma |
| `replace_missing_S_with_clim_S` | True | Fill missing salinity with WOA23 climatology |
| `exclude_high_latitude_profiles_from_clim_cost` | True | Exclude profiles poleward of `dubious_clim_lat_threshold` from climatology cost |
| `dubious_clim_lat_threshold` | 60 | Latitude threshold for the above |
| `clim_cost_qc_with_high_lat_sigma` | True | Step 06: still apply the climatology-cost tests poleward of the threshold to profiles with a high-latitude σ (Arctic file or southern constant) |
| `use_high_lat_sigma` | True | Step 04: use the high-latitude uncertainty files (`high_lat_sigma_file_dict`: `Theta_rmse_se`, `S_rmse_se` from `profile_data/high_latitude_uncertanties/`) for |lat| > `high_lat_sigma_min_lat` (60) |
| `high_lat_sigma_max_dist_km` | 300 | North: use the nearest file cell with a value within this distance, else the CTD σ |
| `southern_sigma_const_dict` | `{'prof_T': 10.0, 'prof_S': 5.0}` | South of −60 (no southern file yet): constant σ at all depths (°C, psu); `None` = CTD σ and no cost tests, as before |
| `remove_extraneous_depth_levels` | False | Step 08: True = trim trailing depth levels with no data; False = keep all files on the full 97-level input depth grid (empty levels = fill, zero weight) |
| `distance_tolerance` | 5000 m | Radius for sub-daily co-location test (step 09) |
| `closest_time` | 120000 (noon) | Preferred sampling time for daily decimation (HHMMSS) |
| `method` | 1 | Decimation method (step 09) |

---

## Step 2 — Map app

The map app is **standalone**: it needs only its three files and a data store built from
post-NCEI profile files. Users with someone else's chain output, or a ready store, don't
need to run the chain. See the **"What you need"** section of
`profile_analysis/README_maps.md` (inputs, commands, options) for full documentation.

The short version:

```bash
# Build the data store (once, or after re-running the NCEI chain).
# --source_root is the run folder; its yearly/ sub-folder is used automatically.
python profile_analysis/build_allsource_bin_timeseries_with_anomalies.py \
    --source_root /path/to/profile_files_NCEI_processed_<YYYYMMDD_HHMMSS>

# Run the app
python profile_analysis/app_profile_map.py
# Open http://127.0.0.1:8050
```

The earlier globe apps (`app_allsource_globe.py`, `app_allsource_anomaly_globe.py`,
`app_profile_scatter.py`) are kept at git tag `globe-apps-2026-10-19`.

---

## Python dependencies

```
numpy
xarray
netCDF4
scipy
matplotlib
cartopy
pandas
pyarrow
dash
plotly
```

The PDF versions of the diagnostics reports use Google Chrome (headless).

---

## Authors

Sweet Zhang, Ian Fenty, Bruce Laughlin  
ECCO Group, JPL / MIT
