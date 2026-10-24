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
| KI-7 | `steps/step06.py` | Dates with day `00` (`YYYYMM00`; year and month known) fail the date test: **18,374 CTD profiles, 1.6% of CTD, all 1992–2002** (up to 7.3% in 2002). Option noted: assign mid-month (the 15th). | Yes; they are zero-weighted, as they always have been | decision needed (mid-month option recorded, not implemented) | yes |
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

## Fixed

(none yet)
