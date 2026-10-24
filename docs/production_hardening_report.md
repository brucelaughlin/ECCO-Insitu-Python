# Production hardening of the NCEI chain: project report

This report explains what is being done to make the NCEI processing chain production-ready, why,
and how each change was checked. It is written to be read on its own, without the code open. A
new section is added at the end of each phase.

- Branch: `production_hardening` (started from `devel_v3` at commit 7fa2cb5)
- Short plan and progress checklist: `docs/production_hardening_plan.md`
- Started: 2026-10-23

---

## Summary for a quick briefing

The NCEI chain works: it processed all 7 sources (about 5.9 million profiles) into model-ready
files. Getting there took a lot of iteration, though. Nearly every new batch of input files turned
up something the code didn't expect, and someone had to read the code, find the problem and patch
it. The goal of this project is a chain that a new user can run on a new set of profiles **without
editing code**, that rejects bad input with a clear explanation instead of crashing, and that has
proper documentation (a Sphinx API reference website) that stays up to date automatically.

The work is split into six phases. The first phase builds a **safety net**: automatic tests that
prove the scientific results do not change while the code is reorganized. That is essential,
because the later phases move and rewrite a lot of code, and we must be able to show that the
processed profiles stay exactly the same.

| Phase | What | Status |
|---|---|---|
| 1 | Safety net: test data and automatic tests | **Done** (2026-10-23) |
| 2 | Configuration file instead of editing code | **Done** (2026-10-23) |
| 3 | Input checks and clear error reporting | **Done** (2026-10-23) |
| 4 | Reorganize into an installable Python package | **Done** (2026-10-23) |
| 5 | Docstrings and Sphinx documentation | Next |
| 6 | Automatic checks on every change (CI), changelog, versioning | Planned |

---

## 0. Starting point (audit, 2026-10-23)

Before planning, I measured the state of the code. The chain's core is 34 Python files and about
9,000 lines: the controller `NCEI.py`, steps 01–13, the shared `tools.py`, the climatology
builders and the `diagnostics/` folder.

### What was found

| Area | Finding | Why it matters |
|---|---|---|
| Documentation in the code | 104 of 279 functions have a docstring (37%). Steps 01–10 and `tools.py` have no description at the top of the file. | An API reference is generated from docstrings, so most of it would be empty today. New readers can't tell what a function expects or returns. |
| Hard-coded paths | `NCEI.py` has 8 paths like `/Users/brucel/ecco/yip/...` for the grid, climatology and uncertainty files; `run_ncei.sh` and the diagnostics have more. | Nobody else can run the chain without editing the code. |
| Settings inside the code | All processing options (floors, thresholds, high-latitude σ, gamma, ...) are variables in the middle of `NCEI.py`. | Changing a setting means editing code. There is no record, next to a run's output, of which settings produced it. |
| Debugging leftovers | `import pdb` in `NCEI.py`, step06, step09 and `tools.py`. | Harmless, but a sign of in-progress code; it can pause a run if a breakpoint is left in. |
| Packaging | No list of required packages and versions; the scripts sit in the repository root and find each other through `sys.path.append(...)` with a fixed path. | It can't be installed; another user can't easily reproduce the software environment. |
| Tests | None. | No way to confirm that a change didn't alter the results, except rerunning everything and comparing by hand. |
| Handling bad input | Each input file goes through the 10 steps inside one try/except. A problem with a file shows up as a crash in whichever step first trips over it. | The error message points at the symptom (for example step 10), not the cause (a variable missing from the input). |

### Lessons from the earlier runs

Looking back at the project log, most of the past iteration came from **input files that differed
from what the code assumed**:

- `prof_Sweight` missing from some files made step 10 crash.
- WOD time values like `SS = 60` and `24:00:00` gave invalid dates.
- Temperature that was already potential temperature, but wasn't recognized as such, was converted a
  second time (the "double conversion" fixed on 2026-10-21).
- Some files had unusual depth grids.

That is why Phase 3 (input checks) is the core of "no tinkering", and why the tests built in
Phase 1 cover those exact cases.

---

## Decisions made (2026-10-23)

| Question | Decision |
|---|---|
| Reorganize the code into a proper Python package? | Yes. |
| Where will the documentation live? | GitHub Pages, since the repository is public. |
| May small test files cut from the real inputs be stored in the repository? | Yes. |
| Where to start? | Phase 1 (safety net). |

---

## Phase 1: safety net (done 2026-10-23)

### Purpose

Phases 2–4 will move settings into a file, add input checks, and reorganize every file into a
package. Each of those changes could, by mistake, alter the processed profiles. Phase 1 builds an
automatic check that answers one question in about 20 seconds: **does the chain still produce
exactly the same output as before?**

There are two kinds of tests:

- A **regression test** runs the whole chain on a small set of input files and compares the
  result with a saved copy of the output (the "reference"). If anything differs, the test fails
  and lists what changed. It doesn't check that the science is right; it checks that the code
  still does exactly what it did before.
- **Unit tests** check single functions against known answers, for example that "24:00:00 on
  31 Dec 2010" becomes "00:00:00 on 1 Jan 2011". They are fast (under 2 seconds in total) and pin
  down the behavior of the parts that caused trouble before.

### 1.1 Test input files ("fixtures")

The full inputs are 246 files and 43 GB, far too large to keep in the repository or to rerun for
every change. So I cut **8 small files**, one per source, from the real input files in
`profile_data/Interp_Profiles`. Each keeps its source file's exact format, variables and
attributes, but only a few hundred profiles.

The profiles in each file were chosen to exercise the parts of the chain that matter:

- **A block of consecutive profiles.** Consecutive profiles often include near-duplicates (the same
  cast reported twice), which step 09's duplicate removal must handle.
- **Profiles north of 60°N and south of 60°S.** These use the high-latitude uncertainties (the
  Arctic rmse files, or the southern constant σ) and get step 06's extra climatology tests.
- **Profiles with out-of-range times** (`SS = 60` and so on), so the date/time fix is exercised.

The selection uses a fixed random seed, so running the cutting script again gives the same files.

| Source | Taken from | Profiles in | Profiles out (after the chain) | Why this file |
|---|---|---|---|---|
| CTD | `WOD_WO_1992_CTD_OSD.nc` | 240 | 88 | Many polar profiles; 24 bad times in the full file |
| Argo (PFL) | `ARGO_WO_2015_PFL_D.nc` | 210 | 210 | Large Argo year with both polar regions |
| XBT | `WOD_WO_1992_XBT.nc` | 190 | 184 | Temperature only (no salinity) |
| Glider | `WOD_WO_2013_GLD.nc` | 200 | 39 | Many bad times; many duplicates removed |
| Mooring | `WOD_WO_2025_MRB.nc` | 150 | 150 | Repeated casts at one place |
| MEOP (seals) | `MEOP_WO_2024_CTD.nc` | 160 | 117 | Southern Ocean; extra variables (O2, CHL, Testim) |
| ITP L2 | `ITP2_WO_2023_CTD.nc` | 120 | 42 | Arctic; the ITP subfolder layout |
| ITP L3 | `ITP_WO_2004_CTD.nc` | 100 | 27 | Arctic; second ITP subfolder |
| **Total** | | **1,370** | **857** | 6.8 MB |

The large drops (CTD, glider, ITP) are the chain's normal quality control, mostly step 09's
duplicate removal; they happen in the full runs too.

The script that cuts the files is `tests/make_fixtures.py`, and the files are in
`tests/fixtures/inputs/`.

### 1.2 Reference outputs

I ran the current code (commit 7fa2cb5, the version that produced the 2026-10-21 run) on the 8
test files, doing the same steps `run_ncei.sh` runs by default:

- steps 01–10 (`NCEI.py`);
- step 11, spatial scaling (monthly windows);
- step 13, packing into ECCO-ready files.

Step 12 (compact) is off by default and isn't included.

The output was saved as the reference in `tests/fixtures/reference/` (7.1 MB):

- `yearly/`: the 8 processed files after step 11;
- `ecco_arrays/`: the 8 model-ready files from step 13;
- `documentation/`: step 13's manifest, the `PROFILES_SIZE.h` / `data.profiles` values, and its
  file-count trade-off table.

**Determinism check.** A regression test only works if the chain gives the same answer every time
it runs. I ran it twice, with different numbers of parallel workers (4 and 2), and compared every
variable in every file: they were identical.

### 1.3 The regression test

`tests/regression/test_chain_regression.py` does the following:

1. It runs steps 01–13 on the 8 test files into a temporary folder, through
   `tests/regression/run_chain.py`.
2. It compares the new output with the reference, file by file, through `tests/regression/compare.py`.
   The comparison checks:
   - the same list of output files;
   - in every NetCDF file: the same dimensions, variables, data types and attributes, and exactly
     equal values (no tolerance, so even a change in the last decimal place is caught);
   - the same JSON manifest and text files.
3. It passes only if nothing differs. Otherwise it lists every difference.

Only two things are deliberately ignored, because they change on every run without meaning
anything:

- the `created` timestamp that step 13 writes into each file;
- absolute folder paths in the manifest, which depend on where the test runs.

**Proof that the test catches real changes.** A test that never fails proves nothing, so I made a
tiny deliberate change: the salinity uncertainty floor from 0.005 to 0.0051 psu (2%). The test
failed and named exactly what changed, for example:

```
ecco_arrays/ITP_L3_001_20040101_20041231.nc:prof_Suncertainty: 1001 values differ (max abs diff 0.0001)
ecco_arrays/ITP_L3_001_20040101_20041231.nc:prof_Sweight: 183 values differ (max abs diff 240.5)
ecco_arrays/ITP_L3_001_20040101_20041231.nc:prof_Scost: 183 values differ ...
```

Then I undid the change, and the test passed again.

**When the external data is missing.** The chain needs the llc90 grid, the WOA23 climatology and
the σ files, about 1.5 GB in folders outside the repository. Where those aren't present (another
computer, or GitHub's automatic checks), the regression test is skipped with a message rather than
failing. Phase 2 makes those paths configurable. Phase 6 then decides how GitHub's servers can run
the full test, for example with trimmed copies of the data files.

### 1.4 Unit tests (31 tests)

| What is tested | File | Why |
|---|---|---|
| Date/time rollover: `SS = 60`, `MM = 60`, `23:59:60`, `24:00:00`, 31 Sep, year-end, leap day; valid times unchanged; fill values and impossible dates left alone | `tests/unit/test_date_time.py` | WOD rounding produced invalid times; this fix was added on 2026-10-14. |
| How step 05 reads the temperature label (`long_name` / `standard_name`): "potential", "in situ", "in-situ", `sea_water_temperature`, no label | `tests/unit/test_step05_T_label.py` | This decides whether T is converted; getting it wrong caused the double conversion. |
| Step 05 does not convert T labeled potential, and obeys a forced "don't convert" setting; the decision is written into the file | `tests/unit/test_step05_T_label.py` | Guards the double-conversion fix directly. |
| The adiabatic temperature gradient matches the published UNESCO 1983 check value (S = 40, T = 40 °C, P = 10,000 dbar → 3.255976 × 10⁻⁴ °C/dbar) | `tests/unit/test_step05_potential_T.py` | Checks the potential-temperature formula against an independent published number. This is a check of correctness, not just of consistency. |
| Step 11's monthly and 10-day windows, including the first and last days of a year and a leap year | `tests/unit/test_step11_windows.py` | Window boundaries decide which profiles count together in the spatial scaling. |

More unit tests will be added in Phase 3 alongside the new input checks.

### 1.5 How to use the tests

From the repository folder, in the `ecco_env1` environment:

```
python -m pytest                 # everything (about 20 s)
python -m pytest tests/unit      # only the fast unit tests (about 2 s)
```

A pass looks like `32 passed`. A failure prints what differs.

**When a change to the results is intended** (for example a new σ rule), the reference must be
regenerated on purpose:

```
python tests/regression/make_reference.py
```

The reason must be recorded in the changelog (from Phase 6) and in this report, and the new
reference committed with the change. That way every change to the science output is deliberate and
documented, never accidental.

### 1.6 Limitations of Phase 1

- **Consistency, not correctness.** The regression test proves the output stays the same, not that
  it is right. The correctness checks so far are the UNESCO value and the earlier comparison with
  the original MATLAB chain (see the project log).
- **Not every situation is covered.** 1,370 profiles can't cover everything that 5.9 million can
  (rare input oddities, very large files). Phase 3's input checks address the rare cases; new
  fixtures will be added whenever a new kind of problem turns up.
- **Step 12 (compact) is not tested,** since it is off by default.
- **Skipped where the external data is missing,** as explained above.

### 1.7 Files added in Phase 1

| File | Purpose |
|---|---|
| `tests/make_fixtures.py` | Cuts the 8 test input files from the full inputs |
| `tests/fixtures/inputs/` | The 8 test input files (6.8 MB) |
| `tests/fixtures/reference/` | The reference outputs (7.1 MB) |
| `tests/regression/run_chain.py` | Runs steps 01–13 the way `run_ncei.sh` does |
| `tests/regression/compare.py` | Compares two output folders exactly |
| `tests/regression/make_reference.py` | Regenerates the reference (only for intended changes) |
| `tests/regression/test_chain_regression.py` | The regression test |
| `tests/unit/*.py` | The 31 unit tests |
| `tests/conftest.py`, `pytest.ini` | Test setup |
| `.gitignore` | Lets the test `.nc` files be stored (other `.nc` files stay ignored) |
| `docs/production_hardening_plan.md` | Short plan and progress checklist |
| `docs/production_hardening_report.md` | This report |

No chain code was changed in Phase 1.

**Commit:** `b78995a` on `production_hardening`, pushed to GitHub.

---

## Phase 2: configuration (done 2026-10-23)

### Purpose

Before this phase, running the chain on another computer, with other data, or with a different
setting meant **editing the Python code**:

- 8 file paths in `NCEI.py`;
- the default folders in `run_ncei.sh`;
- the step 13 settings in `step13_pack_ecco_arrays.py`;
- more paths in the diagnostics and the climatology builders.

That has three problems:

- Every user's copy of the code differs.
- A typo in a setting gives no warning.
- The output doesn't record which settings produced it.

Phase 2 moves **every path and every option** into settings files. The code no longer contains any
machine-specific path, and no setting has to be changed in the code.

### 2.1 How the settings work

The settings are YAML files (a plain-text "name: value" format) in a new `config/` folder. They
are read in layers; each later layer overrides the ones before it:

| Order | Where | What it holds | Stored in git? |
|---|---|---|---|
| 1 | `config/defaults.yaml` | **Every** setting, with its default value and a comment explaining it. This is also the reference list of what can be set. | Yes |
| 2 | `config/local.yaml` | This computer's data locations: where the grid, climatology, σ files, inputs and outputs are. Each user writes their own once, starting from `config/local.example.yaml`. | **No** (it's personal) |
| 3 | `$ECCO_INSITU_CONFIG` | Optional: a settings file named by an environment variable, e.g. for a shared server. | — |
| 4 | `--config my_run.yaml` (`run_ncei.sh -C`) | Settings for one run: any subset, e.g. a different σ floor. | Up to the user |
| 5 | `--set key=value` (`run_ncei.sh -s`) | One setting changed from the command line, e.g. `-s chain.new_floor.prof_S=0.006`. | — |

Paths in `local.yaml` can be relative to one folder, `data_root`, so a new user usually only needs
to change that one line if their data is laid out the same way.

The defaults are **exactly the values the code had before**, so a run with no extra settings
behaves exactly as before. That is what the regression test checks (section 2.4).

### 2.2 Protection against mistakes

| Mistake | What happens now |
|---|---|
| Misspelled setting, e.g. `new_flor` | The run stops before processing: `unknown setting "chain.new_flor" (see config/defaults.yaml for the valid ones)`. Before, a misspelled variable in the code was silently ignored. |
| Wrong kind of value, e.g. `apply_gamma_factor: 3` | The run stops: `"chain.apply_gamma_factor" = 3 has the wrong type (expected bool, like the default True)`. |
| A data file or folder is missing or the path is wrong | The run stops **before processing any file** and lists every missing item, e.g. `paths.grid_dir: does not exist: /...`. Before, this showed up as a crash inside a step, once per file. |
| Relative path but no `data_root` | The run stops with a message saying so. |

### 2.3 A record of the settings with every run

Every run now writes the complete merged settings it used to
**`<run folder>/documentation/run_config.yaml`**. The file starts with the list of settings files
that were merged, in order.

Any output folder therefore records exactly how it was made: which climatology, which σ files,
which floors and thresholds. This matters for comparing runs, and for answering "what settings
produced the files we gave the modelers?" months later.

To see the merged settings without running anything:
`python config_loader.py [--config FILE] [--set key=value]`.

### 2.4 How it was checked

1. **The regression test from Phase 1 passes against the original Phase 1 reference outputs,
   unchanged.** I checked the new outputs against the reference files committed at the end of
   Phase 1 (commit b78995a), not against a newly made copy. Every variable in every output file is
   bit-for-bit identical. The only new file is `documentation/run_config.yaml`, which the reference
   now includes. This is the main proof that moving the settings changed no results.
2. **The test now drives the real entry point.** In Phase 1 the test called the three Python steps
   itself. It now runs `run_ncei.sh`, the same script used for production runs, so the shell
   script and the way settings pass through it are tested too.
3. **Settings really reach the chain.** I ran the chain with
   `-s chain.new_floor.prof_S=0.0051` (a 2% change). The test reported 19 differences (the S
   uncertainties, weights and costs, plus the changed `chain` section in `run_config.yaml`), and the
   saved `run_config.yaml` showed `prof_S: 0.0051` with the `--set` listed as a source. So an
   override isn't silently dropped.
4. **13 new unit tests for the settings loader**
   (`tests/unit/test_config_loader.py`):
   - the defaults load;
   - later layers win, and unchanged values are kept;
   - a misspelled key is rejected;
   - wrong types are rejected (four kinds);
   - `convert_insitu_T_to_potential` accepts `auto`/true/false;
   - relative paths are resolved, or rejected when there is no `data_root`;
   - missing data files are all listed;
   - the saved `run_config.yaml` reads back the same;
   - the example `local.yaml` only uses valid settings.
5. **Every script still starts** (`-h` works for all 10 command-line tools). The two run checkers,
   `check_run_files.py` and `check_ecco_arrays.py`, pass on the test run's output. The diagnostics
   find the same default folders as before.

**Total: 45 tests, all passing (about 21 s).**

### 2.5 A bug found on the way

`build_woa23_climatology.py` (which builds the WOA23 climatology) called
`step06.calculate_potential_T`, but that function lives in **step05**. The script would have
crashed the next time anyone rebuilt the climatology. The reference must be left over from an
earlier renumbering of the steps. It is fixed: the script now calls step05's function, which is the
same UNESCO routine the chain applies to the observations, as the script's description says.
The existing climatology file is not affected; it was built before the renumbering.

### 2.6 What changed, file by file

| File | Change |
|---|---|
| `config/defaults.yaml` (new) | Every setting with its default and an explanation |
| `config/local.example.yaml` (new) | Template for a machine's data locations |
| `config/local.yaml` (new, not in git) | This machine's data locations (`data_root: /Users/brucel/ecco/yip`) |
| `config_loader.py` (new) | Reads, merges, checks and saves the settings; `--get` gives single values to shell scripts |
| `NCEI.py` | The 70-line block of hard-coded paths and options is replaced by reading the settings. Checks the data files exist before starting. Saves `run_config.yaml`. New `--config` / `--set` options; `-i` is now optional (default `run.input_dir`). The `sys.path` line with a fixed path is replaced by the script's own folder. |
| `run_ncei.sh` | Defaults come from the settings instead of fixed paths. New `-C` (settings file) and `-s` (one setting) options, passed to every step. The worker count is now passed to step 11 too. |
| `step11_spatial_scaling.py` | Window and worker defaults from the settings; `--config` / `--set` |
| `step13_pack_ecco_arrays.py` | Packing settings from the `pack` section instead of constants in the file; `--config` / `--set` |
| `optimize_profile_files.py` | No fixed grid path; `--grid_dir` defaults to the settings |
| `build_woa23_climatology.py` | Output folder from the settings or `--out_dir`; the step05 fix above |
| `prebake_woa23_climatology.py` | Defaults from the settings |
| `diagnostics/lib/common.py`, `costs_compute.py`, `costs_scaling.py`, `report_spatial_scaling.py`, `check_run_qc.py`, `report_costs.py` | Default folders and files from the settings |
| `README.md`, `diagnostics/README.md` | New "Settings" section; the parameter table now gives setting names; no machine-specific paths |
| `tests/regression/run_chain.py` | Runs `run_ncei.sh`; can pass extra settings |
| `tests/regression/compare.py` | Compares `run_config.yaml` setting by setting, ignoring machine-specific paths |
| `tests/unit/test_config_loader.py` (new) | The 13 loader tests |
| `tests/fixtures/reference/documentation/run_config.yaml` (new) | The settings of the reference run |
| `.gitignore` | `config/local.yaml` is not committed |

**No step's processing code (steps 01–10) was changed.**

### 2.7 For a new user, after Phase 2

1. Copy `config/local.example.yaml` to `config/local.yaml` and set `data_root` (and any paths that
   differ).
2. Run `./run_ncei.sh`. If anything is missing, the message names it.

### 2.8 Not done yet (later phases)

- **Input checks (Phase 3).** The settings and data paths are checked up front, but the
  **input profile files** are not yet. A file missing a variable can still fail inside a step.
- **Interactive tools.** `profile_analysis/` (the map app and its store builder) and
  `geodesic_binning/` still take their paths as command-line options or constants; they aren't part
  of the chain and will be handled with the packaging (Phase 4).
- **Running the full regression test on GitHub's servers** needs copies of the grid and
  climatology files there (Phase 6).

---

## Phase 3: input checks and clear error reporting (done 2026-10-23)

### Purpose

This phase is the core of "no tinkering for new datasets". As section 0 showed, most past
iteration came from **input files that differed from what the code assumed**: a missing variable,
odd times, unlabeled temperature. The chain found out only when a step crashed partway through,
with an error that described the symptom, not the cause.

Phase 3 makes the chain:

1. **Write down what it needs** from its inputs: the *input contract*, `docs/input_contract.md`.
2. **Check every input file against it before processing anything.** That takes about 12 seconds
   for all 246 input files, since only the small variables are read.
3. **Skip a file that can't be processed, with the reason**, and carry on with the rest, instead
   of failing partway through.
4. **Say what happened to every file** in a run summary.

It also addresses the user's question about climatologies: the chain now accepts any climatology
that meets a written contract, checks it, and offers to build the standard one when none is set.

### 3.1 The input contract

`docs/input_contract.md` lists, for profile files:

- **Required:** what must be there for a file to be processed. A missing item is an *error*, and
  the file is skipped:
  - the `iPROF` / `iDEPTH` dimensions;
  - positions, dates and times;
  - a usable depth grid (meters, increasing);
  - at least one of T and S, each with its weight variable.
- **Conventions** the chain assumes: fill value −9999 or NaN, potential temperature, psu.
- **Warnings:** things the chain handles, but a user should know about. The file is processed and
  the affected profiles are zero-weighted or corrected:
  - bad positions;
  - dates before 1950, impossible, or with no day (`YYYYMM00`);
  - times like `24:00:00`;
  - values outside plausible ranges (often a unit problem);
  - empty variables;
  - unlabeled T.

The checks are in a new module, `validate_inputs.py`. They can also be run alone on a folder,
without running the chain: `python validate_inputs.py <input_dir>`.

### 3.2 What the check found in the real inputs

On all 246 files of `Interp_Profiles`, **every file passes**; 120 have warnings:

| Warning | Files |
|---|---|
| `prof_S` exists but has no valid values (T-only data in a file with an S variable) | 77 |
| Times rolled forward (`SS = 60`, `24:00:00`, ...) | 28 |
| **Day missing** (`YYYYMM00`: year and month known, day `00`): **18,377 profiles**, almost all CTD 1992–2002 | 12 |
| No position (`lon = 99999`) | 7 |
| Latitude outside −90–90 | 7 |
| T values outside −3–40 °C (a handful per file, up to 61 °C) | 6 |
| Longitude out of range | 1 |

**The dates with a missing day are a finding worth noting.** Some CTD profiles have a year and month
but **no day**: the day is `00`, as in `20020100` (January 2002). They have no time either
(`HHMMSS` = 0), and their `prof_date` falls on the last day of the previous month (what a day of 0
works out to).

Step 06 zeroes their weights as "bad dates", so **they never reach the model**. That was true before
this project too; the new check just makes it visible.

How many:

- **18,374 of the 1,174,838 CTD profiles (1.6%)**, all in 1992–2002;
- plus 3 XBT profiles (1992) out of about 720,000.

Within those years the share rises steadily, from about 1–3% in the early 1990s to **7.3% in
2002**; none after 2002:

| Year | CTD profiles with no day | CTD profiles that year | Share |
|---|---|---|---|
| 1992 | 1,290 | 61,551 | 2.1% |
| 1993 | 783 | 63,417 | 1.2% |
| 1994 | 1,653 | 56,509 | 2.9% |
| 1995 | 1,993 | 61,684 | 3.2% |
| 1996 | 1,742 | 50,411 | 3.5% |
| 1997 | 1,539 | 52,840 | 2.9% |
| 1998 | 1,476 | 51,720 | 2.9% |
| 1999 | 1,968 | 51,885 | 3.8% |
| 2000 | 1,887 | 47,039 | 4.0% |
| 2001 | 2,173 | 45,461 | 4.8% |
| 2002 | 1,870 | 25,698 | 7.3% |

The sharp stop after 2002 suggests a few contributing datasets that report only the month; WOD's
metadata could identify them.

**Possible change, not implemented (the user's decision, 2026-10-23):** assign mid-month (the 15th,
with a flag). The profiles would then be kept in the right month, and the monthly spatial scaling
and the model's monthly comparison would be barely affected. Step 09's duplicate removal compares
profiles on the same day, so it would need checking. **This is recorded for the team to decide; the
chain does not change these profiles.**

### 3.3 What happens to each file, and the run summary

Each input file now ends in one of four states, recorded in the new
**`<run>/documentation/run_summary.txt`**:

| State | Meaning |
|---|---|
| `WRITTEN` | processed; number of profiles written |
| `SKIPPED` | failed the input check; the reasons are listed |
| `NO DATA LEFT` | processed, but QC removed every profile (not an error), and the step where it happened |
| `FAILED` | an error during processing, with the step number and the error message (the full traceback is in the log) |

Each file's warnings are listed under it. The summary starts with the counts, for example
`files found: 8 / skipped by the input check: 0 / written: 8 / ...`.

**Exit codes**, so scripts and schedulers can tell what happened:

| Code | Meaning |
|---|---|
| 0 | all good |
| 3 | the run finished, but some files failed (see the summary). Steps 11–13 still run on the files that succeeded. |
| 2 | bad options or settings (nothing processed) |

Before, a failed file was just a message in the long log, and the run reported success.

### 3.4 Climatologies: your own file, or the standard one built for you

This answers the question raised during this phase: "is the design flexible in how climatologies
are made?" It now is:

- **A climatology contract** (`docs/input_contract.md`, section 2). The chain needs monthly
  potential temperature and salinity on a regular lon/lat grid, with named variables and
  coordinates. **Any file meeting it can be used**, whether built by our script, made by a user, or
  from another atlas. The file is checked before the run. It is rejected with the reason if, for
  example:
  - the variable names are wrong;
  - it isn't 12 months;
  - the grid is irregular;
  - it is labeled in-situ temperature;
  - the values are in kelvin;
  - it has undeclared fill values such as −9999 that would be read as real temperatures.
- **The standard one built for you.** As the user suggested: if `paths.climatology_file` points
  to a file that doesn't exist, the run stops before processing and says:
  ```
  For the climatology: point paths.climatology_file at your own file, or build the standard
  WOA23 one with:  python setup_climatology.py --build
  ```
  The new `setup_climatology.py`:
  - checks the current climatology and its pre-baked copies (`python setup_climatology.py`);
  - downloads and builds the standard WOA23 climatology exactly as this project did (`--build`):
    monthly fields to 1500 m stitched onto seasonal fields to 5500 m, and in-situ T converted to
    potential T with the chain's own routine;
  - makes the pre-baked copies (`--prebake`).

  I chose an explicit command over building automatically inside a run, because the build
  downloads 1.4 GB. A run shouldn't start a large download without being asked.
- **A safeguard against mixing climatologies.** The pre-baked files (the climatology
  pre-interpolated to each depth grid, for speed) are used *instead of* the climatology file. Before
  this phase, the chain used every pre-baked file in the folder regardless of which climatology it
  came from. So **switching to a new climatology would have silently kept using the old one's
  pre-baked values** for most files. Each pre-baked file records its source, and the chain now
  uses only those made from the current climatology, listing any it ignores.

### 3.5 Smaller robustness fixes

- **Debugging leftovers removed:** `import pdb` in four files, and two `pdb.set_trace()` calls in
  `tools.py` (inside disabled code, but they could have paused a run if re-enabled).
- **Clearer failure messages:** a crash in a step now logs the error type and the full traceback,
  not just the message.
- **No input files pass the check:** the run says so, instead of silently doing nothing.

### 3.6 How it was checked

1. **The regression test passes; science output unchanged.** The step 01–13 outputs of the 8
   fixture files are bit-identical to the Phase 1 reference. The only new reference file is
   `run_summary.txt`.
2. **One unit test per contract rule** (`tests/unit/test_validate_inputs.py`, 27 tests). Each takes
   a valid fixture file and breaks one thing:
   - drops a required variable (5 variables);
   - drops a weight;
   - drops all observed variables;
   - renames a dimension;
   - reverses or negates the depths;
   - adds bad positions and dates;
   - converts T to kelvin;
   - removes the T label;
   - feeds a file that isn't NetCDF.

   It checks the right error or warning is raised. Also checked: a T-only file is accepted, and all
   real fixtures pass. For climatologies: synthetic files with wrong names, 4 months, an irregular
   grid, in-situ labels, kelvin, and undeclared and declared fill values; plus the real WOA23 file.
3. **A broken file in a real run** (`tests/regression/test_bad_input_handling.py`). The 8 fixtures
   plus a ninth file with `prof_lat` removed go through `run_ncei.sh`. The test checks that:
   - the broken file is skipped;
   - the summary names it with the reason (`"prof_lat" is missing`);
   - **the other 8 files' outputs are identical to the reference**, so one bad file affects nothing
     else.
4. **The pre-baked safeguard** (`tests/unit/test_climatology_guard.py`). A pre-baked file whose
   source is a different climatology is ignored, a matching one is used, and an unusable
   climatology stops the run.
5. **The missing-climatology message.** Running with a nonexistent climatology path stopped before
   creating any output and printed the build instructions.
6. **Tools that read the log still work.** The cost report's parser of the chain log
   (`costs_scaling.parse_step_counts`) reads the new log correctly.

**Total: 76 tests, all passing (about 75 s, two of them full chain runs).**

### 3.7 What changed, file by file

| File | Change |
|---|---|
| `validate_inputs.py` (new) | The profile and climatology contracts as checks; command line for checking a folder or a climatology |
| `setup_climatology.py` (new) | Check / build / pre-bake the climatology |
| `docs/input_contract.md` (new) | The written contract, and what the check found in the real inputs |
| `NCEI.py` | Checks the climatology and pre-baked files, then every input file, before processing; skips files that fail; workers report written / no data left / failed with the step and error; writes `run_summary.txt`; exit code 3 when files failed; logs tracebacks; `import pdb` removed |
| `run_ncei.sh` | Exit code 3 from steps 01–10 lets steps 11–13 run and is passed on at the end; exit codes documented |
| `build_woa23_climatology.py` | Its main work is now a `build()` function, so `setup_climatology.py` can call it |
| `config_loader.py` | The missing-climatology message points to `setup_climatology.py` |
| `step06.py`, `step09.py`, `tools.py` | `import pdb` / `pdb.set_trace()` removed |
| `README.md` | New sections: input checks and the run summary; the climatology |
| `tests/unit/test_validate_inputs.py`, `test_climatology_guard.py`, `tests/regression/test_bad_input_handling.py` (new) | The tests above |
| `tests/fixtures/reference/documentation/run_summary.txt` (new) | The reference run's summary |

**No step's processing (steps 01–13) was changed**, so all the science stays as it was.

### 3.8 Follow-up: making the climatology options obvious

At the user's request, the two climatology options ("use your own file" or "build the standard
one with `setup_climatology.py --build`") are now stated where a new user looks first:

- **The settings files.** A comment box directly above `paths.climatology_file` in
  `config/defaults.yaml`, `config/local.example.yaml` and `config/local.yaml` explains:
  - that a climatology doesn't have to exist already;
  - what a user-supplied file must contain;
  - the build, pre-bake and check commands.
- **The README.** A new **Quick start** section near the top. Its step 2 is "Get a climatology:
  you don't need one already", with both options and the commands.
- **The API documentation (Phase 5).** Its front page and quick start will carry the same text.
  This is recorded in the plan as a requirement.

### 3.9 Not done in this phase

- **`print` → `logging`.** This is moved to Phase 4. The cost report reads the chain's log text
  (`step: N`, `valid prof_T profile count ...`), so the log format has to stay the same. And
  Phase 4 rewrites every module's imports anyway, so doing it there touches each file once.
- **Dates with a missing day** (section 3.2): flagged, not changed (mid-month approximation noted, not implemented). Whether to keep those profiles is a
  science decision for the team.
- **Checks inside steps 01–10** (e.g. step 04's σ file covering the grid) still rely on each
  step's own code. The up-front check covers the input files, which is where the past problems were.

---

## Phase 4: an installable package with one command (done 2026-10-23)

### Purpose

Until now the chain was about 40 loose scripts in the repository's top folder. They found each
other through `sys.path` lines, one of them hard-coded to `/Users/brucel/...`. That causes three
problems:

- **It can't be installed.** Another user, or another machine, has to copy the folder and run
  scripts from inside it.
- **It's hard to document.** Sphinx builds the API reference from an importable package, so loose
  scripts give a thin, disorganized reference.
- **Many entry points.** Each script had its own name and options; there was no single place to
  see what can be run.

Phase 4 turns the code into a Python package, **`ecco_insitu`**, installable with `pip`, with a
single **`ecco-insitu`** command. The science code inside is unchanged.

### 4.1 The new layout

| Before | After |
|---|---|
| `NCEI.py` | `ecco_insitu/pipeline.py` |
| `step01.py` … `step10.py`, `step11_spatial_scaling.py`, `step12_…`, `step13_…` | `ecco_insitu/steps/` (same file names) |
| `optimize_profile_files.py` | `ecco_insitu/packing.py` |
| `config_loader.py`, `config/defaults.yaml` | `ecco_insitu/config.py`, `ecco_insitu/defaults.yaml` (the defaults now ship inside the package) |
| `validate_inputs.py` | `ecco_insitu/validate.py` |
| `setup_climatology.py`, `build_woa23_climatology.py`, `prebake_woa23_climatology.py` | `ecco_insitu/climatology/` (`setup.py`, `build_woa23.py`, `prebake.py`) |
| `tools.py`, `clim_interp.py`, `high_latitude_sigma.py`, `csv_to_nc.py` | `ecco_insitu/` (same names) |
| `diagnostics/*.py` and `diagnostics/lib/*.py` | `ecco_insitu/diagnostics/` |
| (none) | `pyproject.toml`: package definition, dependencies, the `ecco-insitu` command |
| (none) | `ecco_insitu/cli.py`: the `ecco-insitu` command |
| (none) | `ecco_insitu/log.py`: logging set-up (section 4.4) |

Files were moved with `git mv`, so each file's history follows it.

`config/local.yaml` (this machine's paths) stays where it was. For a copy installed without the
repository, its location is `~/.config/ecco_insitu/local.yaml`.

### 4.2 Installing and the `ecco-insitu` command

```
pip install -e .            # in the repository folder (once)
ecco-insitu --help          # lists every command
```

| Command | Does |
|---|---|
| `ecco-insitu run` | the whole chain, steps 01–13 (same as `./run_ncei.sh`, same options) |
| `ecco-insitu steps` | steps 01–10 only (what `python NCEI.py` did) |
| `ecco-insitu scale` / `compact` / `pack` | steps 11 / 12 / 13 |
| `ecco-insitu optimize` | the step 13 packing program with all its options |
| `ecco-insitu validate` | check input files, or a climatology, against the contract |
| `ecco-insitu climatology` | check / `--build` / `--prebake` the climatology |
| `ecco-insitu config` | print the merged settings |
| `ecco-insitu check-files`, `check-arrays`, `check-qc` | the run checks |
| `ecco-insitu report-counts`, `report-costs`, `report-scaling` | the reports |
| `ecco-insitu csv-to-nc` | convert WOD CSV exports to NetCDF |

Every command takes `--help`, and those that use settings take `--config` / `--set`. The commands
work from **any folder**, not only from inside the repository.

**The old commands keep working.** `python NCEI.py ...`, `python step11_spatial_scaling.py ...`,
`python diagnostics/check_run_files.py ...` and the others are now three-line wrappers that call
the matching `ecco-insitu` command. Anyone's existing scripts and habits don't break, and they can
switch when convenient. `run_ncei.sh` works as before and now runs the package's modules.

### 4.3 Dependencies are written down

`pyproject.toml` lists what the chain needs: numpy, xarray, netCDF4, scipy, pandas, pyyaml, pyproj,
geopy, pymatreader and matplotlib. It also defines optional extras:

- `diagnostics`: cartopy, pyarrow, for the maps and reports;
- `test`: pytest;
- `docs`: sphinx.

`pip install -e ".[diagnostics,test]"` installs them. Before, a new user had to discover the
packages by running scripts until one failed to import. The lower version bounds are the major
versions known to work, from the environment the chain was validated in. They are not tested
minimums; Phase 6's automatic tests will show whether they hold.

### 4.4 Logging

The controller's own messages now go through Python's standard `logging`:

- the climatology and pre-baked checks;
- the input-file check;
- the run summary;
- the start-up lines.

Messages that need attention (an ignored pre-baked file, a climatology warning) are logged as
*warnings*. A program using the package can now redirect, filter or silence them with one line, and
by default they print exactly as before.

**Deliberately unchanged:** the per-file output of steps 01–10 (`step: N`,
`valid prof_T profile count ...`, the `!!!` step 05 boxes). Each worker captures it per file and
prints it as one block. The cost report (`report-costs`) reads those exact lines to rebuild how many
profiles each step removed, so their format must not change. Converting the roughly 250 `print`
calls inside the steps to logging would risk breaking that for little benefit. It can be done step
by step later, with the regression test guarding the output.

### 4.5 How it was checked

1. **The regression test passes against the Phase 1 reference.** All step 01–13 outputs of the 8
   fixture files are bit-identical. The reference files were not regenerated in this phase, so this
   is a direct comparison with the outputs of the code before any hardening work.
2. **The new command, run from another folder.** `ecco-insitu run` was started from `/tmp` (not the
   repository). Its outputs were compared with the reference: identical.
3. **The old entry point.** `python NCEI.py` (now a wrapper) ran steps 01–10 on the fixtures and
   wrote the same files.
4. **Every command starts** (`tests/unit/test_cli.py`, 34 tests):
   - all 15 `ecco-insitu` subcommands show their help;
   - the top-level help lists them all, plus the climatology hint;
   - an unknown command gives exit code 2;
   - **all 17 old script names** still run, each started from an empty temporary folder, so nothing
     depends on the current directory.
5. **The log is unchanged where it matters.** The log of a run before and after Phase 4 has the
   same structure: 80 `step:` lines, 80 count lines, 8 file headers, 4 controller lines. The cost
   report's log parser extracts **identical** per-step profile counts from both.
6. **The internal helper modules** (`costs_compute`, `counts_index`, `counts_maps`,
   `counts_timeseries`) also start with `python -m ecco_insitu.diagnostics.<name>`.

**Total: 111 tests, all passing (about 85 s).**

### 4.6 What changed, file by file

- **Moved** (history kept): the 39 files in the table in section 4.1.
- **Imports:** in every moved file, rewritten from `import tools` to
  `import ecco_insitu.tools as tools` and so on. This was done by a script, so every one was caught.
  The code that uses them is unchanged.
  - All `sys.path` manipulation inside the package was removed, including the hard-coded one in
    `profile_analysis/plot_bin_occupancy_per_month.py`, which now finds the package relative to its
    own location.
- **Command-line entry points:** each tool's `if __name__ == '__main__':` block became a `cli(argv)`
  function (the same code), so `ecco-insitu` can call it.
- **New:** `pyproject.toml`, `ecco_insitu/__init__.py` (package description and version
  `1.0.0.dev0`), `ecco_insitu/cli.py`, `ecco_insitu/log.py`, a description file for each
  sub-package (`steps/`, `climatology/`, `diagnostics/`), and the 17 compatibility wrappers.
- **`ecco_insitu/config.py`:** defaults read from inside the package; the location of `local.yaml`
  is decided as in section 4.1.
- **`ecco_insitu/pipeline.py`:** controller messages through logging (4.4).
- **`run_ncei.sh`:** runs `python -m ecco_insitu....` modules; works without installing, by
  pointing Python at the checkout.
- **Docs:** `README.md` (new Repository structure and Installing sections; Quick start uses
  `ecco-insitu`), `diagnostics/README.md` (commands), `docs/input_contract.md`, the `defaults.yaml`
  comments, `docs/profile_file_optimization.md`: all now give the new command names.
- **Tests:** imports updated to the package names; `pytest.ini` folded into `pyproject.toml`; new
  `tests/unit/test_cli.py`.
- **`.gitignore`:** the build folders that `pip install` creates.

**No step's processing code was changed.**

### 4.7 Not done in this phase

- **The interactive map app and its store builder** (`profile_analysis/`) and the old
  `geodesic_binning/` tools are not part of the package. They are separate tools with their own
  dependencies (Dash, Plotly, zarr). They could become an optional extra later.
- **`print` → logging inside the steps** (see 4.4).
- **Publishing to PyPI** (so others could `pip install ecco-insitu` without the repository) isn't
  planned yet. Installing from the GitHub repository works:
  `pip install git+https://github.com/brucelaughlin/ECCO-Insitu-Python`.

---

## Glossary

| Term | Meaning |
|---|---|
| Regression test | An automatic check that the program still gives exactly the same output as a saved reference. |
| Unit test | An automatic check of one function against a known answer. |
| Fixture | A small input file kept for testing. |
| Reference output | The saved output of the code as it was before the changes; the regression test compares against it. |
| Deterministic | Gives exactly the same answer every time it runs on the same input. |
| CI (continuous integration) | GitHub runs the tests automatically on every change and marks the change red or green. |
| Sphinx | The standard tool for Python documentation websites; it builds the API reference from the docstrings. |
| Docstring | The description written at the top of a function: what it does, its inputs and outputs. |
| Package | Code organized so it can be installed with `pip install` and imported anywhere, instead of loose scripts. |
| Input contract | The written list of what an input file must contain for the chain to process it. |
| Pre-baked climatology | A copy of the climatology already interpolated onto one input depth grid, so step 03 doesn't have to do it per file. |
| Exit code | A number a program returns when it ends; scripts use it to tell success (0) from failure. |
| Editable install (`pip install -e .`) | Installs the package as a link to the repository folder, so code changes take effect without reinstalling. |
| Logging | Python's standard way for a program to report messages, which the program that uses it can redirect or filter. |
