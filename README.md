# ECCO In-Situ Processing Pipeline

Authors: Sweet Zhang, Ian Fenty, Bruce Laughlin  
Originally: Sweet Zhang, 2023-05-21  
Updated: 2026-10

---

## Overview

This repository contains a Python pipeline for post-processing in-situ ocean
profile data (from WOD, Argo, MEOP, mooring, and other sources) into a form
suitable for assimilation into the ECCO ocean state estimate. The pipeline:

1. Converts raw WOD CSV exports to NetCDF (`ecco-insitu csv-to-nc`)
2. Runs a 10-step quality-control and preprocessing chain on each profile file (`ecco-insitu steps`)
3. Applies spatial density scaling to the weights (step 11)
4. Packs the result into model-ready profile files that minimize the model's
   memory (step 13, `ecco-insitu pack`; plain-language explanation in
   [`docs/profile_file_optimization.md`](docs/profile_file_optimization.md))
5. Checks and reports on the output (`ecco-insitu check-*` / `report-*`, see [`docs/diagnostics.md`](docs/diagnostics.md))
6. Aggregates the processed profiles into geodesic bins for visualization, with an
   interactive globe app (`profile_analysis/`)

A support data package (LLC90 grid, WOA23 climatology, sigma files, geodesic
bin files) is required and is available on the project Google Drive:
https://drive.google.com/drive/folders/17h0qMS7vVimet8FXieGP1mWhnqnY0ljr

## Quick start

0. **Install** (once, in the repository folder): `pip install -e .`
1. **Point the chain at your data.**
   ```bash
   cp config/local.example.yaml config/local.yaml
   ```
   Then set `data_root` (and any paths that differ) in `config/local.yaml`.
2. **Get a climatology: you don't need one already.** The chain compares every profile with a
   monthly climatology (`paths.climatology_file`). Either:
   - **use your own:** point `paths.climatology_file` at it. Any NetCDF file with
     `potential_T_monthly` and `S_monthly` on (12 months, depth, lat, lon) and a regular lon/lat
     grid works (rules: [`docs/input_contract.md`](docs/input_contract.md), section 2); or
   - **build the standard one:** this project's WOA23 climatology, made on the fly:
     ```bash
     ecco-insitu climatology --build      # downloads WOA23 (~1.4 GB), builds it where paths.climatology_file points
     ecco-insitu climatology --prebake    # optional: faster step 03
     ecco-insitu climatology              # check the setup
     ```
   If the climatology is missing or unusable, a run stops before processing and says so.
3. **Check your input files** (optional; the chain does this anyway):
   ```bash
   ecco-insitu validate /path/to/inputs
   ```
4. **Run:**
   ```bash
   ecco-insitu run
   ```
   The output, the settings used (`documentation/run_config.yaml`) and a per-file summary
   (`documentation/run_summary.txt`) go in a new timestamped folder.

---

## Repository structure

The chain is the Python package `ecco_insitu` (in `src/`), run with the `ecco-insitu` command:

```
ECCO-Insitu-Python/
  pyproject.toml                  — package definition, dependencies, the ecco-insitu command
  config/
    local.example.yaml            — template for this machine's data locations (copy to local.yaml)
  src/ecco_insitu/                — the package
    cli.py                        — the ecco-insitu command (one subcommand per tool)
    run.py                        — ecco-insitu run: the whole chain, steps 01–13, with its log
    pipeline.py                   — steps 01–10 on a folder of input files (ecco-insitu steps)
    steps/                        — step01.py … step10.py, step11_spatial_scaling.py,
                                    step12_profiles_compact.py, step13_pack_ecco_arrays.py
    defaults.yaml                 — every setting, with defaults and explanations
    config.py                     — reads, checks and saves the settings
    validate.py                   — the input-file and climatology contracts (docs/input_contract.md)
    climatology/                  — setup.py (check / --build / --prebake), build_woa23.py, prebake.py
    packing.py                    — the step 13 packing program: any LLC grid / tile size; trade-off report
    clim_interp.py                — climatology interpolation (step 03)
    high_latitude_sigma.py        — high-latitude T/S uncertainties (step 04)
    tools.py                      — shared utilities
    csv_to_nc.py                  — convert WOD CSV exports to NetCDF
    log.py                        — logging set-up
    diagnostics/                  — checks and reports on chain output (docs/diagnostics.md)
  tests/                          — regression and unit tests (python -m pytest)
  docs/                           — input contract, diagnostics, file-packing guide, hardening plan and report
  patches/                        — saved, tested changes not yet applied (see HANDOFF.md)
  profile_analysis/               — the interactive map app and its data store (not part of the package)
    app_profile_map.py, map_common.py, build_allsource_bin_timeseries_with_anomalies.py, README_maps.md
```

### Installing

```bash
pip install -e .                 # in the repository folder; -e keeps it linked to this checkout
pip install -e ".[diagnostics,test]"   # also the map/report extras and pytest
ecco-insitu --help               # all commands
```

---

## Step 0 — Preprocessing: `ecco-insitu csv-to-nc`

Converts a WOD `.csv` export to a set of per-year NetCDF files, with basic
validation of flags, coordinates, and units. Profiles failing validation are
excluded and logged.

```bash
ecco-insitu csv-to-nc -i <input_csv_dir> -d <output_dir>
```

Output filenames: `[stem]_[year].nc`

---

## Step 1 — NCEI processing chain: `ecco-insitu run`

Processes each input NetCDF file through 10 sequential steps (`ecco-insitu steps` runs only these):

| Step | Function | Description |
|------|----------|-------------|
| 01 | `update_prof_and_tile_points` | Interpolates profiles onto LLC90 grid; assigns tile/grid indices |
| 02 | `update_spatial_bin_index` | Assigns geodesic bin IDs (`prof_bin_id_a` at 10242-bin, `prof_bin_id_b` at 2562-bin) |
| 03 | `update_monthly_mean_TS_clim` | Interpolates WOA23 T/S climatology to each profile's location, depth, and day-of-year |
| 04 | `update_sigmaTS` | Assigns T/S uncertainty σ (weight = 1/σ²) from the CTD sigma files; north of 60°N from the high-latitude uncertainty files (`high_latitude_sigma.py`), south of 60°S a constant σ_T = 10 °C, σ_S = 5 psu; records `prof_Vuncertainty_source` (0 CTD, 1 Arctic file, 2 southern constant) |
| 05 | `update_prof_insitu_T_to_potential_T` | Converts in-situ to potential temperature when needed. Decided per file: the `prof_T` metadata (long_name / standard_name) decides if present (all of `Interp_Profiles` is labelled potential); with no label, a check against the climatology decides (below 1000 m in-situ T sits above the potential climatology, growing with depth); if that is unclear too, T is converted. The check never overrides a label (it warns instead). Label-less decisions and disagreements are flagged in the log in a `!!!` box, and every output file records the decision in `prof_T:T_conversion` |
| 06 | `update_zero_weight_points` | Zeros weights on profiles failing quality criteria; optionally excludes high-latitude profiles from climatology cost. Runs before γ, so its cost tests use (ΔV/σ)² |
| 07 | `update_remove_zero_weighted_profiles` | Removes profiles with all-zero T and S weights |
| 08 | `update_remove_extraneous_depth_levels` | Removes trailing depth levels with no valid data (off by default: `remove_extraneous_depth_levels = False`) |
| 09 | `update_decimate_subdaily_profiles` | Decimates sub-daily sampling at the same location to once-daily |
| 10 | `update_gamma_factor` | Multiplies T/S weights by the grid-cell area factor γ (`prof_area_gamma`); moved after QC (was step 05) on 2026-10-12 |

### Post-chain steps (run by `ecco-insitu run` across the whole run directory)

| Step | Command alone | In `ecco-insitu run` | Description |
|------|--------|---------|-------------|
| 11 | `ecco-insitu scale` | on (`--no-scaling` to skip) | Spatial density scaling of T and S weights |
| — | `ecco-insitu report-scaling` | off (`--plots` to run) | Step 11 diagnostic figures |
| 12 | `ecco-insitu compact` | off (`--compact` to run) | Rewrites files compactly in place |
| 13 | `ecco-insitu pack` | on (`--no-pack` to skip) | Packs the yearly files into ECCO-ready profile files (below) |

Defaults come from the `run.*` settings.

**Output layout** (inside each timestamped run folder):

```
profile_files_NCEI_processed_<YYYYMMDD_HHMMSS>/
  yearly/          one file per source and year (steps 01–10, 11, 12)
  ecco_arrays/     step 13: the model-ready profile files only
  documentation/   run_config.yaml (the settings used), run_summary.txt (what happened to each input
                   file), and the step 13 summary: ecco_arrays_manifest.json / .csv, optimize_sweep.csv,
                   PROFILES_SIZE_and_data.profiles.txt (NOBSGLOB / NFILESPROFMAX, data.profiles list)
```

The log of everything printed goes to `z_logs/ncei_log_<YYYYMMDD_HHMMSS>.log` (setting `run.log_dir`).

```
```

**Step 11 — spatial density scaling.** Downweights densely sampled regions so
each (geodesic bin, window) contributes roughly one profile's worth of forcing.
For each variable V in `SCALED_VARS` (`T`, `S`), step 11 counts the profiles
with at least one nonzero `prof_Vweight` in each 10242-bin cell per window
(calendar month by default; `-w 10day` or `spatial_scaling.window: 10day` for 10-day windows), across all sources
combined, and multiplies each profile's `prof_Vweight` by 1/n_V. T and S are
counted separately, so T-only profiles (e.g. Samoa moorings) do not reduce the
S weights of other profiles in the same bin. Files carrying only some of the
weights are scaled for what they have. Adds `prof_Tspatial_scaling_factor` /
`prof_Sspatial_scaling_factor` and renames files to `__ncei_step_11.nc`.
Refuses to run on files that are already scaled — run it on a step-10 copy.
Rewriting the files is nearly all of its run time, so that part runs in parallel
(8 workers by default; `-n` to change).

The original MATLAB (`update_weights_based_on_spatial_density_TRYTOIMPLEMENT.m`, in the git history)
counts all profiles together and applies one factor to both T and S (a Python
version of that, `step11_spatial_scaling_combined_LEGACY.py`, is in the git history
up to commit 1aa5bd1). Unlike the MATLAB, step 11 does not zero weights outside a
fixed range of scaling years.

```bash
ecco-insitu scale <run_dir> --dry_run       # preview
ecco-insitu scale <run_dir>                 # apply in place
ecco-insitu report-scaling <run_dir>        # 28 figures (14 per T/S)
```

Figures go to `<processed_root>/figures/figures_<YYYYMMDD_HHMMSS>/`, outside
the run directory.

**Step 13 — ECCO-ready profile files** (`ecco-insitu pack` with the chain's `pack` settings;
the program itself is `ecco-insitu optimize`, `src/ecco_insitu/packing.py`). New to this? Start with
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
ecco-insitu optimize <profiles_dir>                         # report only
ecco-insitu optimize <profiles_dir> --write <out_dir>       # write the files
ecco-insitu optimize <dir> --nobsglob 35000 --year_rule none
ecco-insitu optimize <dir> --llc 270 --tile 30 --grid_dir <llc270 grid dir>
```
Check with `ecco-insitu check-arrays <run_dir>`.

### Checking a run

Two quick checks should pass on every new run before its output is used (details
in [`docs/diagnostics.md`](docs/diagnostics.md)); both print PASS / FAIL and exit 0 on PASS:

- **`ecco-insitu check-files`: are the yearly files complete and well-formed?**
  Every input has an output or a logged failure; all files at the same step and on
  one depth grid; weights finite, ≥ 0 and never on fill values; every profile used;
  bins and dates valid; and the step 11 factors sum to exactly 1 in every bin-month,
  for T and S. About 15 s.
- **`ecco-insitu check-arrays`: are the model-ready files a faithful repacking?**
  Every yearly profile appears exactly once with identical data; one provenance
  per file, in time order, with whole-year labels; no tile over `NOBSGLOB`; the
  `PROFILES_SIZE.h` values in `documentation/` match the files. About a minute per GB.

```bash
RUN=<run.output_base>_<ts>          # the run folder
ecco-insitu check-files $RUN --input_dir <run.input_dir> --log z_logs/ncei_log_<ts>.log
ecco-insitu check-arrays $RUN
```

A slower scientific check (`ecco-insitu check-qc`, ~1 h) and the reports (profile
counts, costs, step 11 figures) are described in [`docs/diagnostics.md`](docs/diagnostics.md).

### Settings: `config/`

All paths and options live in YAML files, not in the code:

| File | What it is |
|------|-----------|
| `src/ecco_insitu/defaults.yaml` | Every setting, with its default and a comment explaining it. Don't edit it for one machine or one run. |
| `config/local.yaml` | This machine's data locations. Not committed: copy `config/local.example.yaml` and set `data_root` (relative paths are taken from it). |
| `--config my_run.yaml` | Settings for one run: any subset of the keys. |
| `--set key=value` | Change one setting, e.g. `--set chain.new_floor.prof_S=0.006`. |

They are merged in that order (later wins), with `$ECCO_INSITU_CONFIG` (a file path) between
`local.yaml` and `--config`. An unknown key or a value of the wrong type stops the run with a message
naming it, so a typo can't be silently ignored. Each run saves the merged settings it used to
`<run>/documentation/run_config.yaml`.

To see the merged settings: `ecco-insitu config [--config FILE] [--set key=value]`.

### Input checks and the run summary

Before processing anything, the chain checks:

- the settings and data paths;
- the climatology;
- **every input file**, against the input contract in `docs/input_contract.md`: required
  variables and dimensions, depth grid, positions, dates, value ranges, and the T label.

A file that can't be processed is **skipped** with the reason, and the rest of the run continues.
Each run writes `<run>/documentation/run_summary.txt`: per input file, whether it was written,
skipped (and why), lost all its data in QC, or failed (at which step, with the error), plus any
warnings.

`ecco-insitu run` exit codes:

- `0`: all good;
- `3`: the run finished, but some files failed (see the summary);
- `2`: bad options or settings.

To check a folder of inputs without running the chain: `ecco-insitu validate <input_dir>`.

### Climatology

Step 03 needs a monthly potential-T / S climatology: `paths.climatology_file`. Use your own file
(any file meeting the contract in `docs/input_contract.md`, section 2) or build the standard WOA23
one:

```bash
ecco-insitu climatology            # check the climatology in the settings and say what to do
ecco-insitu climatology --build    # download WOA23 (~1.4 GB) and build it where paths.climatology_file points
ecco-insitu climatology --prebake  # optional: pre-baked copies per input depth grid (faster step 03)
```

The chain uses only the pre-baked files that were made from the current climatology.

### Required support files

Set in `config/local.yaml` (`paths` section):

| Setting | Description | Location in support package |
|----------|-------------|---------------------------|
| `paths.grid_dir` | LLC90 grid directory | `grid_llc90/` |
| `paths.sphere_bin_dir` | Geodesic bin assignment files | `grid_llc90/sphere_point_distribution/` |
| `paths.climatology_file` | WOA23 full-depth T/S climatology (NetCDF) | `woa23_climatology/woa23_decav91C0_TS_clim_potential_T_1deg_fulldepth.nc` |
| `paths.prebaked_clim_dir` | Folder of pre-baked climatologies (optional, faster) | `woa23_climatology/` |
| `paths.sigma_files` | CTD T/S uncertainty fields (binary) | `CTD_sigma_TS/` |
| `paths.high_lat_sigma_files` | Arctic T/S uncertainty files (when `high_latitude_sigma.enabled`) | `high_latitude_uncertanties/` |

The chain checks these exist before processing anything, and lists any that are missing.

Pre-baked climatology files (optional but strongly recommended for performance)
are picked up automatically: every `*_prebaked.nc` in `paths.prebaked_clim_dir`
is loaded and matched to input files by the depth grid
stored in the file (`obs_depth`), so adding a new grid needs no code change.
Generate them once with:

```bash
python -m ecco_insitu.climatology.prebake        # defaults: run.input_dir and paths.climatology_file
python -m ecco_insitu.climatology.prebake --input_dir <dir> --source <woa23_fulldepth.nc> --out_dir <dir>
```

This pre-interpolates the WOA23 climatology onto each distinct observation depth
grid found in the input files, eliminating per-profile vertical interpolation at
runtime (step 03). Without pre-baked files the chain still runs correctly but
is significantly slower.

### Running the chain

`ecco-insitu run` timestamps the output folder and writes a log of everything:

```bash
ecco-insitu run                              # input run.input_dir, output run.output_base_<timestamp>
ecco-insitu run --config my_run.yaml         # settings for this run from a file
ecco-insitu run --set chain.new_floor.prof_S=0.006   # change one setting for this run
ecco-insitu run -i /path/to/input            # another input folder
ecco-insitu run -n 8                         # 8 parallel workers
ecco-insitu run -i /path/to/input -n 8 -d /path/to/output   # an exact output folder (must not exist)
ecco-insitu run --no-scaling                 # skip step 11 spatial scaling
ecco-insitu run -w 10day                     # step 11 with 10-day windows
ecco-insitu run --plots                      # also make the step 11 figures
ecco-insitu run --compact                    # also run step 12 compaction
ecco-insitu run --no-pack                    # skip step 13 ECCO-ready packing
```

Parallelism default: `min(8, cpu_count - 1)` workers. Sequential mode: `-n 1`.
On a modern workstation, ~420 files complete in approximately 10–15 minutes
with 8 workers.

Steps 01–10 alone (no steps 11–13, no log file):

```bash
ecco-insitu steps -i <input_dir> -d <output_dir> [-n <n_workers>] [--config FILE] [--set key=value]
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

Every setting, with its default and an explanation, is in `src/ecco_insitu/defaults.yaml`. The main ones:

| Setting | Default | Description |
|-----------|---------|-------------|
| `chain.llcN` | 90 | LLC grid resolution (90 or 270) for step 01's grid / tile assignment |
| `chain.wet_or_all` | 1 | Step 01 grid placement: 1 = wet points only, 0 = all points |
| `chain.respect_existing_zero_weights` | false | Keep weights that are already zero in the input |
| `chain.new_floor` | T 0, S 0.005 | Minimum σ per variable |
| `chain.apply_gamma_factor` | true | Step 10: multiply weights by the LLC90 area factor γ |
| `chain.replace_missing_S_with_clim_S` | true | Fill missing salinity with WOA23 climatology (only for the potential-T calculation; prof_S is not changed) |
| `chain.convert_insitu_T_to_potential` | `auto` | Step 05: `auto` = metadata decides, then the climatology check, then convert (see step 05 above); `true` / `false` force it. Converting an already-potential T cools it a second time (≈0.4 °C at 4800 m; the bug fixed 2026-10-21) |
| `chain.exclude_high_latitude_profiles_from_clim_cost` | true | Exclude profiles poleward of `dubious_clim_lat_threshold` from the climatology-cost tests |
| `chain.dubious_clim_lat_threshold` | 60 | Latitude threshold for the above |
| `chain.clim_cost_qc_with_high_lat_sigma` | true | Step 06: still apply the climatology-cost tests poleward of the threshold to profiles with a high-latitude σ (Arctic file or southern constant) |
| `high_latitude_sigma.enabled` | true | Step 04: use the high-latitude σ for \|lat\| > `min_lat` (60): north, the Arctic files (`Theta_rmse_se`, `S_rmse_se`); south, `south_constant` |
| `high_latitude_sigma.max_dist_km` | 300 | North: use the nearest file cell with a value within this distance, else the CTD σ |
| `high_latitude_sigma.south_constant` | T 10.0, S 5.0 | South of −60 (no southern file yet): constant σ at all depths (°C, psu); null = CTD σ and no cost tests |
| `chain.remove_extraneous_depth_levels` | false | Step 08: true = trim trailing depth levels with no data; false = keep all files on the full 97-level input depth grid |
| `chain.distance_tolerance` | 5000 m | Step 09: radius for the sub-daily co-location test |
| `chain.closest_time` | 120000 (noon) | Step 09: preferred sampling time for daily decimation (HHMMSS) |
| `chain.method` | 1 | Step 09: decimation method |
| `spatial_scaling.window` | monthly | Step 11: profile-density windows (monthly or 10day) |
| `pack.*` | v4r6 values | Step 13: grid, tile size, year rule, NOBSGLOB choice |

-----------|---------|-------------|
| `llcN` | 90 | LLC grid resolution (90 or 270) for step 01's grid / tile assignment |
| `wet_or_all` | 1 | Grid interpolation: 0 = nearest wet point, 1 = all points |
| `respect_existing_zero_weights` | False | Preserve pre-existing zero weights |
| `new_floor_dict` | T=0, S=0.005 | Minimum weight floors |
| `apply_gamma_factor` | True | Apply area-based gamma correction to sigma |
| `replace_missing_S_with_clim_S` | True | Fill missing salinity with WOA23 climatology (only for the potential-T calculation; prof_S is not changed) |
| `convert_insitu_T_to_potential` | `'auto'` | Step 05: `'auto'` = metadata decides, then the climatology check, then convert (see step 05 above); `True` / `False` force it. Converting an already-potential T cools it a second time (≈0.4 °C at 4800 m; the bug fixed 2026-10-21) |
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
