# Model comparison tools: plan and status (work in progress)

Three tasks from the team, 2026-10-23, for checking the processed profiles against the ECCO
model. The code is in `src/ecco_insitu/model_comparison/`. It works and was checked against ECCO
v4r6, but it's still being developed: interfaces may change, and there are no regression tests
for it yet.

| # | Task | Command | Status |
|---|---|---|---|
| 1 | Model output (MDS) → NetCDF | `ecco-insitu mds-to-nc` | **working**, verified bit-identical to the MATLAB |
| 2 | Monthly-mean model fields → `prof_Vestim` | (planned) | **not started**, waiting for monthly-mean fields |
| 3 | Offline cost, exactly as the model | `ecco-insitu model-cost`, `model-cost-bars` | **working**, reproduces the v4r6 bars; final check against a model log pending |

The order chosen was 1 → 3 → 2: (1) was quick, and (3) can be checked against the model straight
away, which then tells whether (2)'s output is right.

## 1. MDS → NetCDF (`ecco-insitu mds-to-nc`)

**The task:** a legacy MATLAB routine turns the model's raw profile output back into NetCDF. It was
"janky, slow, not parallel", one of the last legacy pieces. It was taken to mean the **profiles
package's model equivalents** (`<name>.<tile>.001.equi.data`, one file per tile), not gridded model
diagnostics. That reading was confirmed with the user. The legacy route was
`reconstruct_profile_nc_from_model_data_binary.m`, which calls `MITprof_gcm2nc.m` (MITgcm/MITprof
repository, `MITprof_IO/`). An early draft from another assistant
(`reference/drafts/xmitgcm_gridded_mds_to_nc_draft.py`) used `xmitgcm`, which reads *gridded* MDS
fields, so it doesn't apply to these files.

**Binary layout** (from `MITprof_gcm2nc.m`):

- big-endian float64;
- per profile, a block of (nr+1) × 2 × nvar values, where nr is the number of depth levels: the
  model values, then a mask, then the 1-based profile index at position nr;
- variables in the order T, S, U, V, ptr, ssh, OXY, those present in the input file.

**What the Python version does differently:**

- keeps every variable of the input file (the MATLAB dropped some, including `prof_bin_id_a`);
- checks that each profile appears at most once, and reports profiles with no model record;
- runs datasets in parallel.

**Verified:** all 21 v4r6 iter47 datasets (5.9 M profiles), in about 3.5 minutes. `prof_Testim` and
`prof_Sestim` are **identical** to the MATLAB `_model.nc` files. Six datasets have a few profiles
with no model record (180 in all).

## 2. Monthly means → `prof_Vestim` (planned)

**The task:** from the monthly-mean fields of an ECCO solution (from a run first; PO.DAAC as a next
step), in one command, put the model values at each profile into `prof_Vestim`. It's for cost
checks, so time precision isn't critical.

**Plan:**

- **Horizontal:** the chain uses nearest neighbor (all `prof_interp_weights` are 1), so each
  profile is one LLC90 cell. Profiles map to cells of the compact LLC90 grid by their stored
  `prof_interp_lon` / `prof_interp_lat`; this matched every profile tested.
- **Vertical:** linear interpolation from model levels to the profile depths, following the model's
  `profiles_interp.F`, which has been read.
- **Time:** the profile's month (optionally interpolating between months). The model itself samples
  daily averages, so this will differ somewhat from a run's `.equi` values. That's fine for cost
  checks, and will be stated in the output metadata.
- **Check:** compare with the iter47 `prof_Testim`; close, but not identical.

**Waiting for:** monthly-mean T/S fields from a run (MDS `.data`/`.meta`).

## 3. Offline cost (`ecco-insitu model-cost`)

**The task:** from profile files with weights, `prof_V`, `prof_Vestim` and the geodesic bin label,
compute the cost(s) **exactly as the model does**. The reference is ECCO v4r6's `cost_profiles.F`
(ECCO-GROUP/ECCO-v4-Configurations, `ECCOv4 Release 6/code/`; copy in `reference/mitgcm_cost_code/`,
verified identical to the GitHub version), with v4r6's `PROFILES_OPTIONS.h` (`SAMPLESPLIT_COST`,
`CLIMMASK`) and `data.profiles` (`profiles_mean_indsamples` = 9, all multipliers 1).

**What it computes**, per variable:

- **anomaly cost:** w · mask · (m − d − m̄ + d̄)², per file;
- **mean cost:** w̄ · (m̄ − d̄)², per (geodesic bin, depth), with
  w̄ = 1 / max(s/n, s/9), where s is the bin mean of 1/w.

Here m is the model value, d the data, and the bars are bin means. Details taken from the Fortran:

- the weight is zeroed where `prof_Vclim` ≤ −990 (CLIMMASK);
- the model mask is m ≠ 0;
- depths are matched across files by exact value;
- each bin's mean term counts once.

Output is printed in the model's own `cost_profiles(file,var)=` format, plus `cost_summary.json`.

**Verified:** on v4r6 iter47 (all 21 files), the cost per datum by dataset type matches the team's
"old" bars. For example, the T mean term is 0.355 and S 0.315, and the anomaly terms per type match
to the precision readable from the plots.

**Still to do:** compare with the model's own numbers (STDOUT `cost_profiles(...)=` lines, or the
`costfunction` file) for exact agreement, when a run's output is available.

**Slicing (added 2026-10-25, at the team's request):** besides the model's cost ("final"),
`model-cost` reports the anomaly cost with weight factors removed:

- `no_sf`: weight / spatial scaling ("decimation") factor;
- `no_gamma_sf`: weight / (γ · sf) = (misfit/σ)².

It writes `cost_breakdown.csv`, one row per variable, file, year, latitude band and depth band, with
counts and the three costs. `ecco-insitu model-cost-versions` plots it by type, year, latitude band
or depth band, for one run or several side by side. It works on this package's files (`sf_T`,
`sf_S`) and on older files such as v4r6's (one `prof_spatial_scaling_factor`).

First use: the v4r6 vs 20261021 comparison in
`/Users/brucel/ecco/yip/cost_comparisons/v4r6_vs_20261021_weights/FINDINGS.md`.

`ecco-insitu model-cost-bars` draws the team's bar plots from one or more `cost_summary.json`
files (old vs new).

## To do before these leave "work in progress"

- [ ] Task 2 (monthly means).
- [ ] Task 3: check against a model run's own cost output.
- [ ] Regression tests with small fixtures (a few tiles of `.equi.data`; a small cost case).
- [ ] Settings in `defaults.yaml` where useful (e.g. `profiles_mean_indsamples`).
