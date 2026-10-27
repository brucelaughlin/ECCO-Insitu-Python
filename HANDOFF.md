# Handoff Context — ECCO-Insitu-Python

**Branch:** `model-cost-slicing` (active work, not merged yet; `main` is the default branch. `devel*`,
`production_hardening` and the old experiments were merged or deleted on 2026-10-23). Repository: the fork `brucelaughlin/ECCO-Insitu-Python` (parent
`ifenty/ECCO-Insitu-Python`).
**Documentation website:** https://brucelaughlin.github.io/ECCO-Insitu-Python/ (rebuilt on every push to
`main`). As one PDF: Actions tab → a *tests* run → Artifacts.
**Production hardening (finished 2026-10-23):** plan and progress in `docs/production_hardening_plan.md`;
the full plain-language write-up, for explaining it to others, is `docs/production_hardening_report.md`.
**Last commit:** see `git log`
**Date:** 2026-10-27 (the "Active work (2026-10-27)" section below is the newest; read it first)

## Active work (2026-10-27): WOD flags in the chain, and the deep XBT values (for every agent)

**Read this first if you're another agent working on this project** (e.g. the profile-analysis
framework in `src/ecco_insitu/diagnostics/profile_analyses/`). The detailed write-up, with all the
numbers and scripts, is **`/Users/brucel/ecco/yip/cost_comparisons/wod_flags/README.md`**: use it for
the user's documents and reports. Everything there was checked against the data; inferences are
labeled.

### 1. WOD quality flags are now applied by the chain (commits 2768878, 56ad9ad)

- **The problem:** UCSD's `Interp_Profiles` files carry no per-value quality flags, so values that the
  World Ocean Database flags as bad had full weight. This was the main cause of the higher XBT cost
  of 20261021 vs v4r6.
- **The source:** WOD's standard-level yearly files,
  `https://www.ncei.noaa.gov/data/oceans/woa/WOD/YEARLY/<XBT|CTD|OSD|MRB|GLD>/STD/<DS>S<year>.gz`.
  - **Proof they're the right flags:** v4r6's `prof_Tflag` is exactly WOD's standard-level flag. It
    matches WOD's 2018 export at 100% of 5.35 M levels, and the current files at 99.4% (WOD updated
    some since).
  - **Not the right source:** WOD's per-year observed-level netCDF (`wod_netcdf/`). Its value flags
    remove almost nothing, because UCSD's files seem to drop those values already.
- **The code:**
  - `src/ecco_insitu/wod_flags.py`, and `ecco-insitu wod-flags [--download | --build]`.
  - **Settings:** `chain.apply_wod_flags: true`, `paths.wod_std_dir` (= `/Users/brucel/ecco/yip/wod_std`,
    1.9 GB of files plus a 43 MB `index/`), and `chain.wod_flags_extend_levels`.
  - **When it runs:** before step 01. It writes `prof_Tflag` / `prof_Sflag` by WOD cast number (from
    `prof_descr`), and step 06 test 2 zeroes those values.
  - **Rule:** per-level flags only, as the old chain; not WOD's whole-cast flag.
  - **Depths between WOD standard levels:** the larger flag of the two levels around.
  - **Dependency:** wodpy.
- **Pending the bosses' approval: KI-19, `chain.wod_flags_extend_levels`** (default 2, the user's
  option B).
  - **The rule:** a value beyond the cast's WOD levels takes the end level's flag only if the last 2
    WOD levels at that end are both flagged. `1` = whenever the end level is flagged; `0` = never.
  - **Size:** at most 0.09% of any source's values.
  - **Documented in:** `docs/known_issues.md`, `README.md`, `defaults.yaml`, the module docstring.
- **Coverage** (every WOD input file):
  - **Casts found:** XBT, CTD/OSD and MRB 100%; GLD 99.66% (misses mostly 2019–24).
  - **Values with a real WOD flag:** XBT T 99.77%, CTD 98%, MRB T 98% / S 90%, GLD 97%. The rest are
    beyond the cast's WOD levels, or WOD has no flags for that variable.
- **Values flagged:** XBT T 3.07%, CTD T/S 1.08% / 1.24%, MRB 1.6–1.7%, GLD 1.1–1.4%. Argo, ITP and
  MEOP are not from WOD and get none.
- **Estimated cost per datum** (the 20261021 files with flagged points zeroed; not a rerun):
  - XBT T 0.222 → **0.181** (v4r6 0.170);
  - CTD T 0.050 → 0.048, CTD S 0.063 → 0.060 (the rest of CTD S is the weights);
  - Argo unchanged.
- **Tests:** 120 pass; the regression reference was regenerated deliberately (WOD-flagged points
  changed).
- **A full run with the flags was started 2026-10-27 01:13:**
  `processed_by_NCEI_profile_data/profile_files_NCEI_processed_20261027_011321/`, log
  `z_logs/ncei_log_20261027_011321.log`. Check that it finished: `ecco-insitu check-files`, then
  `check-qc`.

### 2. The deep XBT "spike" in the global mean T: CORRECTION to the profile-analysis notes

- **What the notes say:** the profile-analysis framework's notes (memory
  `project_profile_analysis_framework.md`, bug 5) and its `_MIN_DEPTH_COUNT = 30` rationale say the 7
  XBT values at 2400 m (mean 13.4 °C) are bad data (depth errors), far from climatology but passing
  QC because σ is large. They also say the WOD flags should fix them. **All three are wrong:**
  - **They're real Mediterranean water.** All 7 are in the Mediterranean: six eastern, from 1993 (33–35° N,
    22–32° E), and one at 42° N, 8° E, from 2008. The deep Mediterranean is about 13 °C to the bottom.
  - **They're close to climatology.** |T − T_clim| = 0.01–0.31 °C against the local WOA23 values of
    12.9–13.7 °C. They only look wrong against an open-ocean value (~2.5 °C).
  - **σ is small, not large:** 0.064–0.087 °C. (misfit/σ)² is 0.04–13.3; they pass step 06 because
    that is under its limits (16 for a whole-profile mean, 100 for one point).
  - **The WOD flags (correctly) leave them alone.**
- **Why the spike:** at 2400 and 2600 m the only XBT values left are 6–7 Mediterranean ones, so the
  "global" mean there is a Mediterranean mean: a sampling effect.
  - **Mediterranean share of weighted XBT values:** 6% at 1500–1800 m, 7% at 2000 m, 15% at 2200 m,
    100% at 2400–2600 m.
  - **Report means** (`a01_T_mean_depth`, run 20261021): 4.8 °C at 2200 m (n = 204), 13.4 °C at 2400
    and 2600 m (n = 7, 6), 3.0 °C at 3000 m (n = 13).
  - Masking depths with few values is a fine plotting choice, but the reason should read "sampling
    (Mediterranean only)", not "bad data".
- **Wider check:** all weighted XBT T values deeper than 1500 m and warmer than 6 °C: 2,452 values in
  1,079 casts.
  - 1,608 are Mediterranean, none more than 3 °C from climatology.
  - Most of the rest are Mediterranean outflow (Gibraltar, off Iberia, 6–8 °C) and the Arabian Sea;
    14 are more than 3 °C off.
  - The WOD flags catch 265 of them (11%).
- **The XBT values that ARE suspect are at 3000 m:**
  - **What they are:** 12 casts from two ships: WOD cruise code "JSVY" (1994, south of Tasmania,
    40–62° S) and "JNQN" (1994–95, south of Japan). Each has its deepest value at exactly 3000 m,
    0.6–2.7 °C warmer than climatology.
  - **In WOD's observed levels** (3 casts checked): 34 levels, ending 1500, 2000, 2500, 3000 m, with T
    to 0.1 °C. One is listed as a T-7 XBT (rated to ~760 m), one as type unknown.
  - **Inferred:** these are coarse standard-depth values in WOD's record, not real XBT measurements at
    those depths.
  - **The WOD flags catch only 1 of the 12.**
  - **Possible fix, a team decision, not done:** ignore XBT values below ~2000 m, deeper than any XBT
    probe reaches.

### 3. Scripts and outputs (all in `/Users/brucel/ecco/yip/cost_comparisons/wod_flags/`)

| file | what |
|---|---|
| `README.md` | the full write-up: use for the user's documents and reports |
| `v4r6_vs_wod2018_csv.py`, `std_flags_to_ours.py` | v4r6's flags = WOD's standard-level flags |
| `coverage_all_inputs.py`, `flag_completeness.py` | coverage of every WOD input file |
| `extend_rule_cost.py`, `extend_levels_counts.py` | KI-19: data removed by settings 0 / 1 / 2 |
| `cost_with_wod_flags.py`, `cost_all_types_with_flags.py` | estimated cost effect |
| `deep_warm_xbt.py` (+ `.txt`) | section 2 |

### 4. Next

- Check the 2026-10-27 run; then the team runs the model on its `ecco_arrays/`, and `model-cost` gives
  the real new costs.
- Update the user's report (`cost_comparisons/report_20261025/REPORT.md`) with sections 1–2, or write a
  new short one.
- **Team decisions:**
  - KI-19 (proposed 2);
  - whether to drop XBT values below ~2000 m;
  - XBT no_flag treatment (earlier finding: keep).
- Merge `model-cost-slicing` to `main` (PR on the fork), after the profile-analysis work is committed.

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
    - **DONE (2026-10-27): the WOD flags chain step** (commit 2768878).
      - **Code:** `src/ecco_insitu/wod_flags.py`, `ecco-insitu wod-flags [--download|--build]`.
      - **Settings:** `chain.apply_wod_flags: true` and `paths.wod_std_dir` (= `wod_std/` in
        `local.yaml`).
      - **Data:** files in `/Users/brucel/ecco/yip/wod_std/` (1.9 GB) and the index in
        `wod_std/index/` (43 MB).
      - **What it adds:** WOD's per-level T/S flags on all WOD sources, before step 01. Step 06 test 2
        applies them.
      - **Coverage:** 100% of XBT, CTD and MRB casts; gliders 99.66% (the misses are 2019–24).
      - **Tests:** 118 pass; the regression reference was regenerated deliberately.
      - **Estimated cost** (no rerun, `cost_comparisons/wod_flags/cost_all_types_with_flags.txt`): XBT T
        0.222 → 0.181 (v4r6 0.170); CTD S 0.063 → 0.060; Argo unchanged (not WOD).
      - **Pending the bosses' approval: KI-19, `chain.wod_flags_extend_levels`** (default 2, chosen
        by the user as option B).
        - **What it decides:** values beyond a cast's WOD levels take the end level's flag only if the
          cast's last 2 levels are both flagged; 1 and 0 are the alternatives.
        - **Size:** at most 0.09% of any source's values.
        - **Where it's documented:** `docs/known_issues.md`, `README.md`, the `wod_flags.py`
          docstring and `defaults.yaml`.
      - **Next:**
        - a full chain run (the user's call; long), then ask the team to run the model on it;
        - update `report_20261025` or write a new short report;
        - the unrelated untracked `diagnostics/profile_analyses/` (another agent's work) was left
          uncommitted.
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
