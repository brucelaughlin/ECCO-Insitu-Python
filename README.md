# ECCO In-Situ Processing Pipeline

Authors: Sweet Zhang, Ian Fenty, Bruce Laughlin  
Originally: Sweet Zhang, 2023-05-21  
Updated: 2026-09

---

## Overview

This repository contains a Python pipeline for post-processing in-situ ocean
profile data (from WOD, Argo, MEOP, mooring, and other sources) into a form
suitable for assimilation into the ECCO ocean state estimate. The pipeline:

1. Converts raw WOD CSV exports to NetCDF (`csv_to_nc.py`)
2. Runs a 10-step quality-control and preprocessing chain on each profile file (`NCEI.py`)
3. Aggregates the processed profiles into geodesic bins for visualisation (`profile_analysis/build_allsource_bin_timeseries.py`)
4. Provides an interactive globe app for exploring the data (`profile_analysis/app_allsource_globe.py`)

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
  step01.py … step10.py           — individual pipeline steps
  step11_spatial_scaling.py       — post-chain spatial density scaling (per T/S)
  step11_spatial_scaling_combined_LEGACY.py — old combined-count version (not in chain)
  plot_spatial_scaling.py         — step 11 diagnostic figures
  step12_profiles_compact.py      — optional post-chain compaction
  tools.py                        — shared utilities
  profile_analysis/
    app_allsource_globe.py        — interactive globe app
    build_allsource_bin_timeseries.py  — builds the pickle for the globe app
    README_globe_app.md           — globe app documentation
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
| 04 | `update_sigmaTS` | Assigns T/S uncertainty fields from pre-computed CTD sigma files |
| 05 | `update_gamma_factor` | Applies area-based gamma correction to sigma |
| 06 | `update_prof_insitu_T_to_potential_T` | Converts in-situ temperature to potential temperature; optionally fills missing S with climatology |
| 07 | `update_zero_weight_points` | Zeros weights on profiles failing quality criteria; optionally excludes high-latitude profiles from climatology cost |
| 08 | `update_remove_zero_weighted_profiles` | Removes profiles with all-zero T and S weights |
| 09 | `update_remove_extraneous_depth_levels` | Removes depth levels with no valid data |
| 10 | `update_decimate_subdaily_profiles` | Decimates sub-daily sampling at the same location to once-daily |

### Post-chain steps (run by `run_ncei.sh` across the whole run directory)

| Step | Script | Default | Description |
|------|--------|---------|-------------|
| 11 | `step11_spatial_scaling.py` | on (`-S` to skip) | Spatial density scaling of T and S weights |
| — | `plot_spatial_scaling.py` | off (`-p` to run) | Step 11 diagnostic figures |
| 12 | `step12_profiles_compact.py` | off (`-c` to run) | Rewrites files compactly in-place |

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

The original MATLAB (`update_weights_based_on_spatial_density_TRYTOIMPLEMENT.m`)
counts all profiles together and applies one factor to both T and S; that
behaviour is preserved in `step11_spatial_scaling_combined_LEGACY.py` (writes
a single `prof_spatial_scaling_factor`; skips files lacking `prof_Sweight`).
Unlike the MATLAB, neither Python version zeroes weights outside a fixed range
of scaling years.

```bash
python step11_spatial_scaling.py <run_dir> --dry_run   # preview
python step11_spatial_scaling.py <run_dir>             # apply in-place
python plot_spatial_scaling.py   <run_dir>             # 28 figures (14 per T/S)
```

Figures go to `<processed_root>/figures/figures_<YYYYMMDD_HHMMSS>/`, outside
the run directory.

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
```

Parallelism default: `min(8, cpu_count - 1)` workers. Sequential mode: `-n 1`.
On a modern workstation, ~420 files complete in approximately 10–15 minutes
with 8 workers.

Alternatively, call `NCEI.py` directly:

```bash
python NCEI.py -i <input_dir> -d <output_dir> [-n <n_workers>]
```

### Output

One NetCDF file per input file, written to the output directory under a
subdirectory named by source (e.g. `CTD_WOD/`, `PFL/`). Each output file
contains the original profile data plus added fields including:

- `prof_bin_id_a`, `prof_bin_id_b` — geodesic bin assignments
- `prof_Tclim`, `prof_Sclim` — WOA23 climatology at each profile
- `prof_Tweight`, `prof_Sweight` — quality weights (spatially scaled after step 11)
- `prof_Tspatial_scaling_factor`, `prof_Sspatial_scaling_factor` — step 11 factors
- `prof_Tcost`, `prof_Scost` — climatology misfit cost (computed in step 07 with pre-scaling weights; not updated by step 11)

### Parameter reference

Parameters set in `NCEI.py` under the `NEED PATHS / PARAMETERS` block:

| Parameter | Default | Description |
|-----------|---------|-------------|
| `llcN` | 90 | LLC grid resolution (90 or 270) |
| `wet_or_all` | 1 | Grid interpolation: 0 = nearest wet point, 1 = all points |
| `respect_existing_zero_weights` | False | Preserve pre-existing zero weights |
| `new_floor_dict` | T=0, S=0.005 | Minimum weight floors |
| `apply_gamma_factor` | True | Apply area-based gamma correction to sigma |
| `replace_missing_S_with_clim_S` | True | Fill missing salinity with WOA23 climatology |
| `exclude_high_latitude_profiles_from_clim_cost` | True | Exclude profiles poleward of `dubious_clim_lat_threshold` from climatology cost |
| `dubious_clim_lat_threshold` | 60 | Latitude threshold for the above |
| `distance_tolerance` | 5000 m | Radius for sub-daily co-location test (step 10) |
| `closest_time` | 120000 (noon) | Preferred sampling time for daily decimation (HHMMSS) |
| `method` | 1 | Decimation method (step 10) |

---

## Step 2 — Globe app

See `profile_analysis/README_globe_app.md` for full documentation.

The short version:

```bash
# Build the pickle (once, or after re-running the NCEI chain)
cd profile_analysis
python build_allsource_bin_timeseries.py \
    --source_root /path/to/profile_files_NCEI_processed \
    --output      /path/to/allsource_bin_timeseries.pkl

# Run the app
python app_allsource_globe.py
# Open http://127.0.0.1:8050
```

---

## Python dependencies

```
numpy
xarray
scipy
matplotlib
cartopy
pandas
dash
plotly
```

---

## Authors

Sweet Zhang, Ian Fenty, Bruce Laughlin  
ECCO Group, JPL / MIT
