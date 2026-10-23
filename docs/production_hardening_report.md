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
| 2 | Configuration file instead of editing code | Next |
| 3 | Input checks and clear error reporting | Planned |
| 4 | Reorganize into an installable Python package | Planned |
| 5 | Docstrings and Sphinx documentation | Planned |
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
