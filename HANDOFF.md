# Handoff Context — ECCO-Insitu-Python

**Branch:** `main` (the only branch; `devel*`, `production_hardening` and the old experiments were
merged or deleted on 2026-10-23). Repository: the fork `brucelaughlin/ECCO-Insitu-Python` (parent
`ifenty/ECCO-Insitu-Python`).
**Documentation website:** https://brucelaughlin.github.io/ECCO-Insitu-Python/ (rebuilt on every push to
`main`). As one PDF: Actions tab → a *tests* run → Artifacts.
**Production hardening (finished 2026-10-23):** plan and progress in `docs/production_hardening_plan.md`;
the full plain-language write-up, for explaining it to others, is `docs/production_hardening_report.md`.
**Last commit:** see `git log`
**Date:** 2026-10-23

## Start here (new session)

```bash
cd /Users/brucel/ecco/yip/ECCO-Insitu-Python && git pull     # on main
pip install -e ".[diagnostics,test,docs]"                    # once per environment (ecco_env1 has it)
python -m pytest                                             # ~90 s; must pass (incl. bit-identical regression)
ecco-insitu --help                                           # every command
```

- **Settings:** `src/ecco_insitu/defaults.yaml` (all, with comments) + `config/local.yaml` (this
  machine's paths; not in git).
- **Whole chain:** `ecco-insitu run`.
- **Before changing science:**
  - read `CHANGELOG.md` and `docs/known_issues.md`;
  - an intended change to the outputs needs `python tests/regression/make_reference.py` and a
    **Science:** changelog entry.
- **Every push** runs the tests on GitHub (Actions tab), and **every push to `main`** republishes the
  docs. Use pull requests within the fork; when opening one, check the **base repository** is
  `brucelaughlin/...`, not `ifenty/...` (GitHub defaults to the parent).

---

## Current state

### Inputs (`/Users/brucel/ecco/yip/profile_data/Interp_Profiles/`)
246 files, all on one 97-level depth grid (2–6000 m). Deliberately moved to
`profile_data/Raw_Profile_Data_Unused_By_Ecco/`: Argo A and R files (Argo is
delayed-mode **D only**), PFL_BGC (95% duplicate of PFL T/S), CTD_WOD `*_BGC` (no
T/S), SOCAT (no T/S), MRB_WOD `36_depths/` and `97_depths/` (overlapping
duplicates), Samoa (own deep grid). Two ITP L2 files (2008, 2018) are empty and
three (2012, 2014, 2017) have bad salinity; the chain rejects all five.

### Chain (`ecco-insitu run`: steps 01–10 in `src/ecco_insitu/pipeline.py`, then post-chain steps)

Since the production-hardening work (branch `production_hardening`), the chain is the package
`src/ecco_insitu/` with the `ecco-insitu` command, and all settings are in
`src/ecco_insitu/defaults.yaml` + `config/local.yaml`. File and setting names below that predate it
(`NCEI.py`, `run_ncei.sh`, `-S`, `remove_extraneous_depth_levels = False`) map to: `pipeline.py`,
`ecco-insitu run`, `--no-scaling`, `chain.remove_extraneous_depth_levels: false`.

| Step | What | Notes |
|---|---|---|
| 01–04 | grid/tile points, geodesic bins, WOA23 climatology, σ | step 01 assigns LLC90 cells and 30×30 tiles |
| 05 | potential T | converts only when needed (2026-10-21 fix, below); checks every file against the climatology |
| 06 | QC (zero weights) | runs **before** γ since 2026-10-12, so its cost tests use (ΔV/σ)²; skips the cost tests poleward of 60° |
| 07 | drop all-zero-weight profiles | |
| 08 | trim empty depth levels | **off** by default (`remove_extraneous_depth_levels = False`) → every file on the 97-level grid |
| 09 | decimate sub-daily repeats | removes ~92% of glider and ~61% of ITP profiles |
| 10 | γ (grid-cell area factor) | moved here from step 05 on 2026-10-12 |
| 11 | spatial scaling, T and S separately | rewrites files in parallel (8 workers, `-n`; ~3× faster); `--no-scaling` to skip |
| 12 | compaction | off by default (`--compact`) |
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

### Temperature double conversion: found and fixed (2026-10-21)
**All runs before 2026-10-21 have T too cold at depth and must not be used.** The raw
files (`Interp_Profiles`) already hold *potential* temperature (`prof_T` long_name says
so); step 05 converted it again. Effect: −0.005 °C at 40 m, −0.04 °C at 400 m, −0.2 °C at
2800 m, −0.4 °C at 4800 m. Salinity was never affected. It showed up as a large
ECCO v4r6 time-mean T cost (good per-profile "wiggles", bad bin means).
- **Fix (`step05.py`, `NCEI.py`):** `convert_insitu_T_to_potential = 'auto'`. The file's
  metadata decides (potential → leave; in-situ → convert). With no label, a check against
  the climatology decides; if that is unclear too, T is converted.
- **The check never overrides a label.** Label-less decisions and disagreements print a
  `!!!` box under `step: 5` in the log. Every output file records the decision in
  `prof_T:T_conversion`.
- **How the check works:** below 1000 m, in-situ T sits above the (potential) WOA23
  climatology by an amount growing with depth. Tested on CTD, Argo and gliders (potential)
  and on synthetic in-situ versions (in-situ); XBT and MEOP 2010 are "unclear" (too few
  deep data).
- **Verified:** run 20261021_100417 — every profile's T equals the raw T in all 7 sources
  (5.93 M profiles, largest difference < 0.0001 °C); `check_run_files.py` and
  `check_ecco_arrays.py` PASS. QC kept 826 more profiles than before (less deep bias).
- **ECCO v4r6 cost check** (`/Users/brucel/ecco/yip/cost_comparisons/ecco_v4r6_old_vs_new_20261021/`):
  `cost_profiles.F` re-implemented offline with the iteration-47 model equivalents
  (`profile_data/from_OU_for_comparison/`). On profiles both pipelines share, our corrected
  T reproduces the old pipeline's deep CTD time-mean cost within 2.7% at every depth; the
  double conversion files were up to 40× higher. Slides (PNG / PDF) are in that folder.
- **Not explained:** XBT T differs from the old pipeline's XBT files by ~0.03–0.1 °C (our T
  matches our raw files; likely a different XBT input version).

Date/time clean-up (in `NCEI.py`, before step 01): out-of-range `HHMMSS` /
`YYYYMMDD` values from WOD rounding (seconds or minutes = 60, 24:00:00, day past
month end) are rolled forward with real date arithmetic, so they agree with
`prof_date` (≈3,900 profiles, mostly gliders).

### Step 13 / `optimize_profile_files.py` — model-ready files
Target: ECCO v4r6 (LLC90, 30×30 tiles; 117 tiles, 113 active). The model
reserves per tile `NFILESPROFMAX` × `NOBSGLOB` profile slots (× levels), where
`NOBSGLOB` is the most profiles any one file puts in any one tile. The program
minimizes that product: one provenance per file (ITP L2/L3 separate, XBT flag +
noflag together), chronological, files labeled in whole calendar years (several
whole years, or parts of one too-dense year with the year's label). Default
`NOBSGLOB = 'auto'` = fewest files within 5% of the minimum product. Needs only
the LLC size and tile size — not blank tiles or bathymetry (the model ignores
profiles on land / in blank tiles itself). Writes the `PROFILES_SIZE.h` values and
the `data.profiles` file list. **Plain-language guide:
`docs/profile_file_optimization.md`.** Validate with
`diagnostics/check_ecco_arrays.py <run>`.

Result for the 2026-10-21 run: 18 files, `NOBSGLOB` 43,355, `NFILESPROFMAX` 18 (fits the
v4r6 build: 100,000 and 35). Argo gets 3 files (1997–2013, 2014–2021, 2022–2026).

### Runs
- **`profile_files_NCEI_processed_20261021_100417`** — **the one to use**: step 05 fix,
  high-latitude σ (Arctic file + southern constant), new layout. 5.94 M profiles, 7
  sources, 241 files. `check_run_files.py` and `check_ecco_arrays.py` PASS
  (`z_logs/check_ecco_arrays_20261021_100417.log`).
- `..._20261014_170404` — the **double conversion files** (T too cold at depth; given to
  the team for the first cost comparison). Kept only for that comparison; do not use.
- All older runs were deleted.

### Diagnostics (`ecco-insitu check-*` / `report-*`, see `docs/diagnostics.md`)
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

0. **New since 2026-10-23 (production hardening, model-comparison tools):**
   - **Model comparison, task 2** (monthly-mean model fields → `prof_Vestim`): waiting for
     monthly-mean T/S fields from the team. Plan in `docs/model_comparison.md`.
   - **Model comparison, task 3:** check `ecco-insitu model-cost` against a model run's own cost
     output (STDOUT `cost_profiles(...)=` lines or the `costfunction` file) when one is available.
   - **Known issues** (`docs/known_issues.md`):
     - *safe fixes, any time:* KI-1 `np.NaN` (csv_to_nc, LLC270 loaders), KI-2 `logical_and` (LLC270),
       KI-10, 11, 13–17;
     - *team decisions:* KI-4 CSV T label, KI-5 CSV years, KI-6 step 06 T ≥ −2 °C, KI-7 the 18,374
       missing-day CTD profiles (mid-month option noted, not implemented), KI-8 (0° N, 0° E),
       KI-9 step 09 HHMMSS distance.
   - **README support-data line:** confirm the Google Drive folder is current and list exactly
     which files the chain needs. The Arctic σ files are probably not in it; the climatology can
     now be built (`ecco-insitu climatology --build`).
   - **Optional:**
     - tag `v1.0.0`;
     - trimmed data so the full-chain regression tests also run on GitHub;
     - a proposal to the team's repository (a PR from the fork's `main` to `ifenty/main`).
   - **Paused:** the deck.gl profile viewer (`docs/deckgl_viewer_plan.md`, with hosting notes).
1. **Hand the 2026-10-21 files to the model team**, and tell them that all earlier runs
   (including the double conversion files) have the cold-T bias.
2. **XBT temperature offset** against the old pipeline (~0.1 °C): compare our raw XBT
   files with theirs (version, fall-rate correction).
3. **Confirm with the model team**: the `PROFILES_SIZE.h` values step 13 prints;
   `NLEVELMAX` could drop from 140 to 97 (our depth grid), shrinking every tile's
   arrays by ~30%.
4. **To-do next time — drop profiles with no geodesic bin** (approved idea, not
   applied). Tested patch: `patches/step06_code12_drop_unbinned_profiles.patch`
   adds QC code 12, toggled by `zero_weights_of_unbinned_profiles` in `NCEI.py`.
   **Note (2026-10-23):** the patch predates the package layout (`src/ecco_insitu/`) and the settings
   file; it no longer applies as is. To use it, re-make it against `src/ecco_insitu/steps/step06.py`,
   with the switch as a `chain.*` setting in `defaults.yaml`.
5. **Southern Ocean uncertainty field**, when the team has one (replaces the constant
   σ south of 60°S, setting `high_latitude_sigma.south_constant`).
6. **Map app** (`profile_analysis/app_profile_map.py`): its store is built from the
   2026-10-21 run (`z_profile_file_analysis/allsource_bin_timeseries_with_anomalies_20261021_100417`).
   The older stores (20260930, 20261014) can be deleted. Rebuild after each chain run with
   `profile_analysis/build_allsource_bin_timeseries_with_anomalies.py --source_root <run>`.
   Possible later feature (paused): model-cost color modes from the `prof_Testim` values.
7. **Does step 13 keep `prof_T:T_conversion`** in the model-ready files? Not checked.
8. **Old map-app stores** (20260930, 20261014) can still be deleted; the old `processed_by_NCEI` runs too,
   if no longer needed.

---

## Key file locations

| What | Path |
|------|------|
| The package | `src/ecco_insitu/` (installed with `pip install -e .`; command `ecco-insitu`) |
| Whole run / steps 01–10 | `src/ecco_insitu/run.py` (`ecco-insitu run`), `src/ecco_insitu/pipeline.py` (`ecco-insitu steps`) |
| Steps | `src/ecco_insitu/steps/step01.py` … `step13_pack_ecco_arrays.py` |
| Settings | `src/ecco_insitu/defaults.yaml` (all, with comments), `config/local.yaml` (this machine) |
| Input checks / climatology | `src/ecco_insitu/validate.py` (`docs/input_contract.md`), `src/ecco_insitu/climatology/` |
| File-packing program | `src/ecco_insitu/packing.py` (`ecco-insitu optimize`; guide: `docs/profile_file_optimization.md`) |
| Checks & reports | `src/ecco_insitu/diagnostics/` (see `docs/diagnostics.md`) |
| Tests | `tests/` (`python -m pytest`) |
| Saved, tested patches | `patches/` |
| Map app (World / Arctic / Southern Ocean) and its store build | `profile_analysis/` (see `README_maps.md`) |
| Input profiles | `/Users/brucel/ecco/yip/profile_data/Interp_Profiles/` |
| Unused inputs | `/Users/brucel/ecco/yip/profile_data/Raw_Profile_Data_Unused_By_Ecco/` |
| Chain output | `/Users/brucel/ecco/yip/processed_by_NCEI_profile_data/` |
| Reports | `/Users/brucel/ecco/yip/cost_reports/`, `/Users/brucel/ecco/yip/insitu_reports/`, `/Users/brucel/ecco/yip/cost_comparisons/` |
| Logs | `z_logs/` |
| WOA23 climatology | `/Users/brucel/ecco/yip/woa23_climatology/` |
| ECCO v4r6 config (tiles, PROFILES_SIZE.h, data.profiles) | github.com/ECCO-GROUP/ECCO-v4-Configurations, `ECCOv4 Release 6/` |

## Architecture notes (map app store)
- Zarr layout: 2-D `(n_months, 97)` arrays per source per bin, with a shared
  `ym_index` (YYYYMM int32).
- Profile records are in SQLite (`_profiles.db`), not the JSON, to keep start-up
  memory manageable.
- The app auto-selects the most recent timestamped store; `--data <path>` overrides.
- `aggregate_bins` uses `ProcessPoolExecutor(fork)`; if it deadlocks, run with
  `OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1`.
