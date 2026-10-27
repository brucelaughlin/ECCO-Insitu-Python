# Known issues

Problems found in the code that are **not fixed yet**. Each has an ID so it can be referred to
in issues, commits and the changelog. When one is fixed, move it to the bottom section with the
version (or commit) that fixed it, and add a line to `CHANGELOG.md`.

Found while documenting the package (2026-10-23; see `docs/production_hardening_report.md`,
section 5.4). **None of these affects current production runs on `Interp_Profiles` with LLC90**
unless stated.

- **Status:** *open*, or *decision needed* when the team has to choose the behavior.
- **Science:** whether a fix would change the processed profiles. If so, the regression reference
  is regenerated on purpose when it is fixed.

## Would crash if the code path were used

| ID | Where | Problem | Affects current runs? | Status |
|---|---|---|---|---|
| KI-1 | `csv_to_nc.py`; `tools.py` LLC270 loaders (`load_llc270_grid_step2/4/5`, lines ~737–840) | Uses `np.NaN`, which numpy 2 removed (installed: 2.4). `ecco-insitu csv-to-nc` and any LLC270 run would crash. | No (CSV conversion and LLC270 aren't used) | open; safe fix: `np.nan` |
| KI-2 | `tools.py` LLC270 loaders (~737, 783, 839) | `np.logical_and(lat == 0, lon == 0, bathy == 0)` passes the third array as numpy's *output* argument, so the blank-cell test is only `lat == 0 & lon == 0`, and it overwrites `bathy` with the result. | No (LLC270 only) | open; safe fix: chain two `logical_and` calls |
| KI-3 | `steps/step01.py`, `step02.py`, `step04.py` | Only LLC90 is fully supported. With another `chain.llcN`, step 01 fails with a NameError, step 02 assumes 90, and step 04 always loads the LLC90 grid. | No (LLC90) | open |

## Could change results (decision needed)

| ID | Where | Problem | Affects current runs? | Status | Science |
|---|---|---|---|---|---|
| KI-4 | `csv_to_nc.py` (~189) | WOD CSV temperature is labeled "potential temperature", but WOD CSV T is normally **in-situ**. Step 05 trusts the label and wouldn't convert. | Only files made by `csv-to-nc`; not `Interp_Profiles` | decision needed | yes, for those files |
| KI-5 | `csv_to_nc.py` (~347, 435, 791) | `unique_years` is shared across all CSV files, and each file's last year is written after the file loop, so a multi-file conversion may lose years. | Only `csv-to-nc` | open | yes, for those files |
| KI-6 | `steps/step06.py` (~119) | The legal T range is −2 to 40 °C. Seawater near freezing can be about −2.1 °C, so the coldest polar water may be zero-weighted. | Yes, possibly a few polar points | decision needed | yes |
| KI-19 | `wod_flags.py` (`level_flags`); setting `chain.wod_flags_extend_levels` | **DECISION PENDING APPROVAL. WOD flags for values beyond a cast's WOD standard levels.** Our input files can have values deeper than the cast's deepest WOD standard level (or shallower than its shallowest), where WOD gives no flag. The original chain never had such values. These are 0.2% (XBT T) to 7% (MRB S) of values. **Proposed (default, `2`):** flag such a value, with the nearest level's flag, only if the cast's last **2** WOD levels at that end are both flagged, so that a single flagged bottom level doesn't remove everything below it. Alternatives: `1` = whenever the end level is flagged; `0` = never. Extra values flagged with `2` (with `1`): XBT T 7,127 (8,452); CTD T/S 5,628 / 7,499 (8,650 / 12,412); MRB T/S 4,201 / 2,231 (6,786 / 5,481); GLD T/S 17,223 / 26,475 (25,877 / 40,116). At most 0.09% of a source's values (MRB S); the flags themselves remove 1.1–3.1%. Counts: `cost_comparisons/wod_flags/extend_levels_counts.txt`. **To change:** set `chain.wod_flags_extend_levels` in a config file, or `--set chain.wod_flags_extend_levels=1`. | Yes, very slightly | decision needed (proposed: 2) | yes, for these values |
| KI-8 | `steps/step06.py` (~220) | A profile at exactly (0° N, 0° E) is rejected as a missing position. The code's own comment questions this. | Possibly a handful | decision needed | yes |
| KI-9 | `steps/step09.py` | "Closest to `closest_time`" is measured on the HHMMSS number, not in seconds (11:59:59 vs 12:00:00 counts as 4041 apart), so among near-noon profiles the choice can be slightly off. The `method` setting is ignored. | Yes, the choice among same-day duplicates | decision needed | yes (which duplicate is kept) |
| KI-10 | `steps/step01.py` (~266, 272) | Any unexpected error is caught and reported as "you may have no valid data", so a real bug looks like an empty file (shown as "no data left" in `run_summary.txt`). | Only if a bug occurs | open | no |

## Minor

| ID | Where | Problem | Status |
|---|---|---|---|
| KI-11 | `steps/step06.py` (~324, 340) | A printed percentage for the conductivity test uses the wrong variable's original weights (log text only). | open |
| KI-12 | `diagnostics/check_run_qc.py` (~92) | Check 3 (comparison with a baseline run) only pairs files named `__ncei_step_10`, so it can't compare two current runs (`__ncei_step_11`). | open |
| KI-13 | `climatology/prebake.py` (~180), `steps/step12_profiles_compact.py`, `csv_to_nc.py` | Stale text: prebake prints old activation instructions (activation is automatic); step 12 stamps "converted by profiles_compact.py"; misspellings, and `units` written as `unit` in `csv_to_nc`. | open |
| KI-14 | `pipeline.py` (~30) | Imports cartopy and matplotlib without using them, which slows start-up and makes cartopy a hard dependency. | open |
| KI-15 | `steps/step13_pack_ecco_arrays.py` (`cli`); `check_ecco_arrays`, `check_run_files`, `step12` | `ecco-insitu pack` always exits 0; a few commands end with `sys.exit` instead of returning an exit code. | open |
| KI-16 | `climatology/prebake.py`, `diagnostics/common.py` | Read the settings when imported, so a broken `config/local.yaml` stops them importing (other modules report it when run). | open |
| KI-17 | `diagnostics/costs_report.py` (`to_pdf`) | The PDF step uses a macOS-only Chrome path. | open |

## Fixed or decided

| ID | Where | Outcome |
|---|---|---|
| KI-18 | input files (`Interp_Profiles`) | **Fixed 2026-10-27: WOD quality flags added by the chain.** The input files carry no per-value quality flags (UCSD confirmed their files never carried WOD's flags), so values that WOD flags as bad had full weight. This was the main cause of the higher XBT cost of the 20261021 run compared with v4r6. The chain now adds WOD's standard-level per-level flags (as the original MATLAB chain had them) before step 01; see `src/ecco_insitu/wod_flags.py` and `cost_comparisons/wod_flags/` for the check against v4r6's flags (99.4% agreement on XBT 1992–1998). |
| KI-7 | `steps/step06.py` | **Decided 2026-10-25: drop them** (no change needed). The 18,374 CTD profiles (1.6% of CTD, 1992–2002) with day `00` fail the date test and are removed by step 07; none is in the output. They aren't in the old v4r6 dataset (only 2 positional near-matches, neither the same cast), the day was already 00 in UCSD's earlier delivery, and 98% sit on a regular grid at 15–33°N, 110–127°E (South China Sea / Luzon Strait / Taiwan): a month-resolution source in WOD, not a processing bug. The mid-month option was considered and not adopted. |
