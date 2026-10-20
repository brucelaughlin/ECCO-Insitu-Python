# Diagnostics: checks and reports on NCEI chain output

Everything you run is in this folder; `lib/` holds the helper modules.
All commands are run from the repo root. Each tool defaults to the newest
finished run (`profile_files_NCEI_processed_<ts>` with every file at step ≥ 11)
unless you pass a run directory.

| Script | Question it answers | Typical time | Output |
|---|---|---|---|
| `check_run_files.py` | Is the run complete and well-formed? | ~15 s | console; exit code 0 = PASS |
| `check_run_qc.py` | Are the values scientifically plausible? | ~1 h | console |
| `check_ecco_arrays.py` | Are the step 13 ECCO-ready files complete, correctly labeled, within the NOBSGLOB cap, and identical to the yearly files? | ~1 min per GB | console; exit code 0 = PASS |
| `report_profile_counts.py` | How many profiles UCSD delivered vs how many survive processing, per source / year / bin? | ~10 min (scan) + plots | `insitu_reports/report_<run>__<ts>/` (HTML, PNGs) |
| `report_costs.py` | How big are the T/S costs vs climatology, per source / year? | ~1 min + report | `cost_reports/cost_<run>__<ts>/` (HTML + PDF, per-profile costs) |
| `report_spatial_scaling.py` | What did step 11 do (counts, occupancy, mean scaling factor, per T/S)? | minutes | `processed_by_NCEI_profile_data/figures/figures_<ts>/` (PNGs) |

Report output folders are under `/Users/brucel/ecco/yip/`. Run folders hold
`yearly/`, `ecco_arrays/` (model-ready files only) and `documentation/` (step 13
summary files); every tool accepts the run folder and reads the per-year files from
`yearly/` (older runs, with source folders at the top level, still work).
`check_ecco_arrays.py` reads the manifest from `documentation/` and flags any
non-profile file left in `ecco_arrays/`.

## After every chain run

```bash
RUN=/Users/brucel/ecco/yip/processed_by_NCEI_profile_data/profile_files_NCEI_processed_<ts>

python diagnostics/check_run_files.py $RUN \
    --input_dir /Users/brucel/ecco/yip/profile_data/Interp_Profiles \
    --log z_logs/ncei_log_<ts>.log

python diagnostics/check_ecco_arrays.py $RUN          # the step 13 model-ready files

nohup python -u diagnostics/check_run_qc.py $RUN --summary > z_logs/qc_<ts>.log 2>&1 &
```

## Each tool

### `check_run_files.py` — structural validation
Checks:
- **Completeness:** all files at the same final step; every input has an output or a logged failure (`--input_dir`, `--log`).
- **File structure:** required variables present; `prof_interp_*` shapes as the Fortran expects; every file on one depth grid (`--allow_mixed_grids` if step 08 trimming is on).
- **Weights:**
  - finite and ≥ 0;
  - no weighted fill (weight > 0 only where the value and its climatology are real);
  - every profile has some nonzero weight.
- **Bins and dates:** bins in range (unbinned land-cell profiles are a warning); valid dates that agree with `prof_date`.
- **Step 11:** scaling factors sum to exactly 1 in every bin-month, for T and S.

Run it on every new run, before using the output. Exit code 0 = PASS; warnings
(e.g. the few hundred profiles on land cells with no bin) don't fail the check.

### `check_ecco_arrays.py` — the step 13 model-ready files
Compares `ecco_arrays/` with the `yearly/` files they were packed from, using the
manifest in `documentation/`. Checks:
- **Nothing lost or duplicated:** every yearly profile appears in exactly one
  packed file, and every per-profile variable is identical to the yearly file.
- **Packing rules:** one provenance per file; profiles in time order; one depth grid.
- **Year labels:** whole-year files never mix part of a year with other years, and
  their labels are whole calendar years (`YYYY0101`–`YYYY1231`).
- **Model limits:** no file puts more than `NOBSGLOB` profiles in any model tile
  (same grid, tile size and positions as step 13 used); the manifest's `NOBSGLOB` /
  `NFILESPROFMAX` match the files.
- **Folder contents:** `ecco_arrays/` holds only profile files.

Run it after every run that includes step 13 (it reads every packed file, so it
takes about a minute per GB). Exit code 0 = PASS.

```bash
python diagnostics/check_ecco_arrays.py <run_dir>
python diagnostics/check_ecco_arrays.py --yearly <dir> --arrays <dir> [--docs <dir>]
```

### `check_run_qc.py` — scientific QC
Checks:
- **Cost outliers** (T and S) and climatology coverage. Profiles the chain's QC exempts from its cost
  tests (poleward of 60° on the CTD σ, i.e. without a high-latitude σ) are not flagged.
- **Positions:** interpolation displacement, interpolated position vs the declared one.
- **Physical consistency:** T/S water-mass envelope; T vs S weight mismatch.
- **Climatology:** WOA23 spot checks at fixed anchor points.
- **Regression:** comparison against an earlier run with `--baseline <old_run>`.

Use `--summary` for compact output. It prints only at the end, so run it with
`nohup ... &` as above.

### `report_profile_counts.py` — UCSD vs processed profile counts
UCSD-website-style (https://ecco.ucsd.edu/assim.html) profiles-per-year charts
per source, kept-fraction heatmaps, a summary table, and bin-count maps (10242
and 2562 bins; all sources and per source; all years and per year).

```bash
python diagnostics/report_profile_counts.py                         # newest run
python diagnostics/report_profile_counts.py --processed <run_dir>
python diagnostics/report_profile_counts.py --reuse_index <report>  # re-plot without rescanning
python diagnostics/report_profile_counts.py --sources ITP XBT_WOD    # quick test
python diagnostics/report_profile_counts.py --per_year_pngs         # also one PNG per year
```
Definitions:
- **UCSD T / S:** the profile has at least one valid T (S) value.
- **processed T / S:** the profile has at least one nonzero weight, i.e. it is actually used. Climatology-filled S has zero weight and isn't counted.
- **effective:** the sum of the step 11 scaling factors over used profiles, roughly the Argo-equivalent count.
- **Bins:** UCSD profiles are binned with the same nearest-LLC90-cell lookup as step 02.

Caveat: the UCSD side is whatever is in `Interp_Profiles` *now*. If files were
moved out after the run was made, that source can show "kept" > 100%.

To make a PDF of the HTML page: open it in Chrome and print, or see how
`lib/costs_report.py:to_pdf` uses headless Chrome.

### `report_costs.py` — cost vs climatology
Computes `cost = (prof_V − prof_Vclim)² × prof_Vweight` for T and S at every
point that enters the fit (value, climatology and weight non-fill; weight > 0).
The final weight is γ / σ² × sf<sub>V</sub>, so every statistic and figure is
shown for **three versions**, to see what each factor does:
- **a)** nothing removed — the cost the model sees (γ · sf<sub>V</sub> · (ΔV/σ)²);
- **b)** step 11 sf<sub>V</sub> removed — γ · (ΔV/σ)²;
- **c)** γ and sf<sub>V</sub> removed — (ΔV/σ)², the misfit in σ units (100 = 10σ).

The explainer at the top says which version the run's QC step tested (it reads
the step order from the run's log: since 2026-10-12 QC runs before γ, so it
tested c).

Report contents:
- summary per version (median, mean, 99th percentile, max, profiles with a point
  ≥ 100, split within ±60° / north of 60°N / south of 60°S), per-source and per-source-year tables,
  distributions, maps of high-cost profiles, the largest outliers with provenance;
- the scaling factors γ, sf<sub>T</sub>, sf<sub>S</sub> (tables and distributions),
  and where and why sf<sub>T</sub> ≠ sf<sub>S</sub> (T-only profiles by source, map);
- profiles removed by QC and by decimation, per source and year (from the log).

Per-profile data for more plots: `costs/<source>/<file>__cost.nc` (per level,
with γ and sf) and `profiles.parquet` (one row per profile). Cost tables show
whole numbers from 1 up and two decimals below 1 (rounded down); the CSVs keep
full precision.

```bash
python diagnostics/report_costs.py                      # compute (~1 min) + report
python diagnostics/report_costs.py --report_only        # re-make report from newest cost dir
python diagnostics/report_costs.py --log <ncei_log>     # if the run's log is not in z_logs/
```

### `report_spatial_scaling.py` — step 11 diagnostics
For T and S, 28 figures:
- profile counts on a 1° grid and per geodesic bin (global, NH, SH);
- P(≥1 profile per window);
- the mean scaling factor.

Also run automatically by `./run_ncei.sh -p`.

```bash
python diagnostics/report_spatial_scaling.py <run_dir> [--window monthly|10day] [--dpi 150]
```

## `lib/`

| Module | Used by |
|---|---|
| `common.py` | paths, latest-run selection, source/subtype parsing, geodesic bin lookup and map rasters (both reports) |
| `counts_index.py`, `counts_timeseries.py`, `counts_maps.py` | `report_profile_counts.py` (scan → per-profile index; time-series figures; maps) |
| `costs_compute.py`, `costs_report.py` | `report_costs.py` (compute + store costs; figures, HTML, PDF) |

The interactive map app (`app_profile_map.py`) and the Argo bin-occupancy plots
(`plot_bin_occupancy_per_month.py`) are in `profile_analysis/`.
