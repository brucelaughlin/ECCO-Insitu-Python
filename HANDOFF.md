# Handoff Context — ECCO-Insitu-Python

**Branch:** `main` (the only branch; `devel*`, `production_hardening` and the old experiments were
merged or deleted on 2026-10-23). Repository: the fork `brucelaughlin/ECCO-Insitu-Python` (parent
`ifenty/ECCO-Insitu-Python`).
**Documentation website:** https://brucelaughlin.github.io/ECCO-Insitu-Python/ (rebuilt on every push to
`main`). As one PDF: Actions tab → a *tests* run → Artifacts.
**Production hardening (finished 2026-10-23):** plan and progress in `docs/production_hardening_plan.md`;
the full plain-language write-up, for explaining it to others, is `docs/production_hardening_report.md`.
**Last commit:** see `git log`
**Date:** 2026-10-26 (the "Active work" section below is the newest; read it first)

## Active work (2026-10-25 / 26): why the 20261021 profile costs differ from v4r6's, and WOD flags

**Branch in progress:** `model-cost-slicing` (pushed, NOT merged). It adds the factor-removed cost
versions to `ecco-insitu model-cost` (writes `cost_breakdown.csv`: by var, file, year, latitude band,
depth band; versions `final`, `no_sf` = weight/sf, `no_gamma_sf` = weight/(γ·sf) = (misfit/σ)²), plus
`ecco-insitu model-cost-versions` (plots). Tests pass (102). Merge with a PR on the fork when ready.

**The team's question:** the model's anomaly cost for our 20261021 profiles is higher than v4r6's for
XBT, ArgoS and CTDS (Ou's bar plots). The boss asked for the cost per datum without decimation (spatial
scaling) and γ ("slice and dice the profile costs offline").

**Data:**
- **v4r6 ("Old"):** the profile files v4r6 has always used, made by the old MATLAB chain (21 files, 21
  Dec 2023): `/Users/brucel/ecco/yip/profile_data/from_OU_for_comparison/input_insitu/`. Model
  equivalents in `.../model_equivalents/run.v4_rls6.102g.iter47/{mds,netCDF}`.
- **20261021 (ours):** `processed_by_NCEI_profile_data/profile_files_NCEI_processed_20261021_100417/`
  (`ecco_arrays/` = the files the model ran on).
  - **Model equivalents from Ou's run `run.v4_rls6.102i9.iter47`:**
    `/Users/brucel/ecco/yip/ecco_v4r6_model_data/model_equivalents/{mds,netCDF}`. On S3 they sit under a
    folder named `..._20261014_170404`, but the profiles match 20261021 exactly (verified).
  - **v4r6 monthly means:** `ecco_v4r6_model_data/{THETA,SALT}_mon_mean` (MDS); for the planned
    monthly-means → V_estim tool, not built yet.
- **Ou's MATLAB converter:** `cost_comparisons/external_scripts/convert_prof_gcm2nc` (with gcmfaces-master
  and MITprof-master). Our `ecco-insitu mds-to-nc` reproduces its NetCDF bit for bit (18 files).
  - Ou's MITprof_gcm2nc has an extra "mean" mode (reads `T_mean.bin` etc., adds
    `prof_V_mean / Vestim_mean / Vweight_mean`); not ported yet.

**Report** (verified twice, all claims rechecked): `/Users/brucel/ecco/yip/cost_comparisons/report_20261025/REPORT.md`
(+ png/, pdf/). It has a Terms section (two kinds of flag; what "reproducing the bars" means) and a list
of corrections. Analyses and scripts: `cost_comparisons/v4r6_vs_20261021_weights/` (`split_by_validity.py`,
`b_points_quality.py`, `noflag_quality.py`, `xbt_qf_vs_noflag.py`, `missing_day_check.py`, ...). Per-run
cost tables: `cost_comparisons/{v4r6_iter47,r20261021}_cost_breakdown/`.

**Findings (confirmed):**
- **Our offline cost:** agrees with Ou's bars for both data sets, to the precision of the plots. It's a port
  of ECCO v4r6 `cost_profiles.F`; the MITgcm version computes the same cost.
- **XBT:** higher in all three versions, so it's NOT decimation or γ.
  - **Main cause:** ~228k XBT points (on casts in both data sets, shared depths) that v4r6 gave zero
    weight because **WOD's per-value quality flag** marked them bad. They have full weight in ours: median
    |T−T_clim| 2.1 °C vs 0.48 °C, 2.2 cost/datum.
  - **Where the old flags came from:** the old chain read them directly from WOD's CSV export ("F"
    column) in `transform_WOD13prof_to_MITprof.m` (`prof_Tflag = prof_T_wod_flag`). The 2018-05-08 WOD
    export is on disk: `sample_data/ecco-insitu/sweet_gdrive/NODC_20180508/ORIGINAL_CSV/`.
  - **UCSD's files never carried WOD flags:** their earlier files had `prof_Tflag` = all zeros
    (placeholders); the current `Interp_Profiles` files have no flag variable. UCSD confirmed this
    (2026-10-26): "used to include a prof_T/Sflag field, but it was just filled with zeros ... I remove
    the empty fields."
  - **Agreement elsewhere:** on points both runs weight, the data agree (ours 0.05–0.09 °C warmer in the
    top 100 m); cost ~10% higher from ~4% larger weights and different bin means.
- **XBT `no_flag` files** (`*_XBT_noflag.nc`, casts submitted to WOD without originator flags) are 64% of
  our XBT data and mostly good (median |T−T_clim| 0.53 vs 0.52 °C for QF'd).
  - **The bad points are in both kinds** at similar rates (1.9% vs 1.5%), so dropping no_flag files is not
    the fix.
  - **Per-value flags exist for them too:** in v4r6's files the no_flag casts have per-value flags at the
    same rate. So "no_flag" refers to originator flags.
  - **A team member's view:** include no_flag (comparable cost, more coverage).
- **CTD S:** higher because of the weights (decimation/γ). **Argo S:** mostly the ~650k extra casts in our
  run (2021–2026). **ITP/MEOP:** lower thanks to the new high-latitude σ.
- **Missing-day CTDs** (18,374, day 00): not in the old dataset; day already 00 in UCSD's earlier delivery;
  98% at 15–33°N, 110–127°E. **Decision: drop them** (the chain already does); KI-7 closed.
- **`scripps_data` folder** (now `/Users/brucel/ecco/yip/x_old_scripps_data_DONOTUSE`): UCSD's previous
  delivery, saved 2026-04-14 (internal dates 2021–Jan 2026). Its flag variables are all zeros.

**Decided approach (the user, 2026-10-26):** don't ask UCSD to reprocess. **Get the modern WOD flags
ourselves.**
- **Source:** WOD per-year netCDF at `https://www.ncei.noaa.gov/data/oceans/ncei/wod/<year>/wod_<inst>_<year>.nc`
  (`xbt, ctd, mrb, pfl, osd, ...`; e.g. wod_xbt_1995.nc ≈ 0.8 GB, updated 2025-06-27).
- **Variables:** `wod_unique_cast`, `Temperature_WODflag` (per value), `Temperature_WODprofileflag`
  (cast), `z_WODflag`, `Temperature_origflag` (originator), `Temperature_row_size` (ragged arrays).
- **Matching:** our casts carry the WOD cast number in `prof_descr` (`WOD-XBT//<cast>:...`).
- **Plan:**
  1. Download WOD XBT 1992–1998.
  2. Match by cast number.
  3. Map WOD flags onto our 97 standard depths. A design choice: e.g. flag a standard-depth value if a
     neighbouring WOD value used in the interpolation is flagged, or the cast's profile flag is bad.
  4. Write `prof_Tflag`; step 06 applies it unchanged.
  5. Recompute the XBT cost: does it drop toward v4r6's 0.17?
  - **Status: test done 2026-10-26.** Full results: `cost_comparisons/wod_flags/README.md`.
    - **v4r6's `prof_Tflag` is WOD's standard-level flag:** it matches the 2018 CSV export at 100% of
      5.35 M levels.
    - **WOD 2025 per-year netCDF:** it has observed-level flags only, plus a per-cast profile flag.
      - **Observed-level flags:** they hit only 3.8 k of our points. Where WOD flags a value, our file
        usually has none; inferred: UCSD's interpolation drops them.
      - **Profile flag:** it catches 72% of the v4r6-flagged points we keep.
    - **XBT 1992–98 cost/datum:** 0.282 → 0.240 with the 2025 profile flag, 0.242 with v4r6's flags
      copied by cast number, 0.229 with both. v4r6 is 0.188. That closes about half the gap.
    - **Standard-level flags found (later, 2026-10-26):**
      - **Source:** `ncei.noaa.gov/data/oceans/woa/WOD/YEARLY/<SRC>/STD/<SRC>S<year>.gz`, read with
        wodpy.
      - **Agreement:** 99.4% with v4r6's flags.
      - **Effect:** they flag 2.3% of our 1992–98 XBT points; cost 0.282 → 0.237. As before, about half
        the gap.
    - **Next:**
      - download STD files for all years and sources, write `prof_Tflag`/`prof_Sflag` in a chain step;
    - **IN PROGRESS (2026-10-26, evening): the WOD flags chain step.**
      - **The user's request:** apply WOD's flags to every WOD source (XBT, CTD/OSD, MRB, GLD; T and S),
        matching the old chain.
      - **Downloads:** done, 1.9 GB, `/Users/brucel/ecco/yip/wod_std/{xbt,ctd,osd,mrb,gld}/<DS>S<year>.gz`,
        1992–2026.
      - **wodpy:** installed in `ecco_env1`.
      - **Written so far:** `src/ecco_insitu/wod_flags.py` (download, index build, check, `apply()`,
        `cli`). It has NOT been run or tested, and is NOT wired in.
      - **Remaining:**
        1. Add `paths.wod_std_dir` (null) and `chain.apply_wod_flags` (false) to `defaults.yaml`, and
           `wod_std_dir: wod_std` to `config/local.yaml` and `local.example.yaml`.
        2. Register the `wod-flags` command in `cli.py` and add `wodpy` to the dependencies in
           `pyproject.toml`.
        3. In `pipeline.py`, call `wod_flags.apply` before step 01, or as part of step 06 test 2. Run
           `check()` up front when the flags are enabled, and pass `wod_std_dir` through the
           function arguments.
        4. Run `ecco-insitu wod-flags --build` to make `wod_std/index/`.
        5. Check: the `std` mask from `cost_comparisons/wod_flags/std_flags_to_ours.py`; for XBT
           1992–98 the flags should hit about 283 k weighted points, 2.28%.
        6. Unit test of `level_flags` and `cast_numbers`; docs (`input_contract`, `known_issues`,
           sphinx api, CHANGELOG).
        7. A full chain run, then `model-cost`.
      - **Note:** the old chain applied only per-level flags, not WOD's whole-cast flag. Our module
        does the same.
      - the team decides which flags step 06 applies, and for which sources.

**Lessons (the user was embarrassed by overstated claims sent to the team and UCSD):** verify every claim
against the data before it goes in a report or message; label inference vs confirmed; don't summarize
from memory. The user prefers doing processing ourselves over asking UCSD to reprocess.

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
