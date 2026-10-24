# Production hardening + API docs: plan and progress

**Resume here.** This file is the single source of truth for this project. Update the progress log
at the bottom after each working step, and commit it with the work.

The full write-up for readers (what was done, why, and how it was checked, in plain language) is
`docs/production_hardening_report.md`. **Add a section to it at the end of every phase.**

- Branch: `production_hardening` (from `devel_v3` at 7fa2cb5)
- Started: 2026-10-23

## Goal
- The NCEI chain runs on new profile sets without code edits.
- It has a Sphinx API doc, hosted on GitHub Pages (the repo is public).
- Docs, changelog and code changes stay in step, enforced by CI.

## Decisions (by the user, 2026-10-23)
1. Restructure into a package: yes.
2. Repo is public, so docs go to GitHub Pages.
3. Small fixture files cut from the real inputs may be committed.
4. Start with Phase 1.

## Phases
1. **Safety net**
   - Regression fixtures: small files per source, cut from `profile_data/Interp_Profiles`.
   - Reference outputs from the current code; a pytest that runs the chain on the fixtures and compares
     the outputs exactly.
   - Unit tests for the tricky pieces.
2. **Configuration**
   - A YAML config for all paths and options; no hard-coded `/Users/...` paths.
   - The config used is saved in each run's output.
3. **Robustness**
   - Validate input files up front against a documented input contract, with clear per-file reasons.
   - `logging` instead of print.
   - Remove the pdb imports.
   - A run summary.
4. **Package**
   - Code moves into `ecco_insitu/`, with a `pyproject.toml`.
   - An `ecco-insitu` command line; `run_ncei.sh` becomes a thin wrapper.
   - The regression test must still pass unchanged.
5. **Docs**
   - NumPy-style docstrings everywhere.
   - Sphinx (autodoc + napoleon): a user guide plus the API reference.
   - **Requirement from the user:** the docs' front page and Quick start must make it
     abundantly clear that a climatology need not exist already: either "use your own file" or
     "`setup_climatology.py --build` makes the standard one". Reuse the README's Quick start and the
     `defaults.yaml` box.
6. **Keep in step**
   - GitHub Actions: tests, a docs build with warnings as errors, a docstring-coverage check, and
     publishing to Pages.
   - `CHANGELOG.md`, a PR template, and a version number stamped in output files.

## Rules
- Science results must not change, unless a change is deliberate, recorded here and in the changelog,
  and the reference outputs are regenerated on purpose.
- One phase per PR or commit series; this file is updated each time.

## Progress log
- 2026-10-23
  - Created the branch and this plan.
  - Phase 1 started: choosing fixtures.
- 2026-10-23, **Phase 1 done**
  - **Fixtures:** `tests/make_fixtures.py` cuts 8 input files into `tests/fixtures/inputs/` (1,370
    profiles, 6.8 MB).
    - One per source: CTD, PFL, XBT, GLD, MRB, MEOP, ITP L2/L3.
    - Each holds a contiguous block plus profiles beyond 60°N and 60°S and bad-HHMMSS profiles; the
      selection is seeded.
  - **Reference:** `tests/fixtures/reference/` (7.1 MB) is the output of the code at 7fa2cb5:
    `yearly/` + `ecco_arrays/` + `documentation/`.
    - Two runs gave bit-identical outputs, so the chain is deterministic.
  - **Regression test:** `tests/regression/test_chain_regression.py` runs steps 01–13 through
    `run_chain.py`, the same steps `run_ncei.sh` runs with its defaults, and compares everything
    exactly with `compare.py`.
    - Ignored: the `created` timestamp and absolute paths.
    - About 20 s.
    - Checked that it fails: changing the S floor 0.005 → 0.0051 made it fail and name the variables
      that changed.
    - Skipped when the external data (llc90 grid) is missing, e.g. in CI, until Phase 2 makes the
      paths configurable.
  - **Unit tests (31):** date/time rollover, step05's T-label decision and forced modes, the UNESCO
    check value for the adiabatic gradient, and step11's windows.
  - **Commands:** `python -m pytest` runs everything; `python -m pytest tests/unit` runs only the fast
    tests.
  - **To regenerate the reference after an intended change:** `python tests/regression/make_reference.py`.
  - pytest was installed into `ecco_env1`.
- 2026-10-23, **Phase 2 done** (details in the report, Phase 2 section)
  - **Settings:** `config/defaults.yaml` (every setting, with today's values), `config/local.yaml`
    (machine paths, gitignored; the template is `local.example.yaml`), `--config` / `--set`
    (`run_ncei.sh -C` / `-s`).
    - `config_loader.py` merges and validates them (unknown key or wrong type → error), checks the data
      paths, and saves `<run>/documentation/run_config.yaml`.
  - **Code changes:** `NCEI.py`, `run_ncei.sh`, steps 11 and 13, optimize, both climatology builders,
    and the diagnostics all read the settings.
    - No hard-coded `/Users` paths remain in the chain or the diagnostics.
  - **Bug fixed:** `build_woa23_climatology.py` called `step06.calculate_potential_T`, which is now in
    step05.
  - **Tests:**
    - The regression test now runs `run_ncei.sh` itself.
    - The outputs are bit-identical to the Phase 1 reference (checked against the committed files);
      only `run_config.yaml` was added to the reference.
    - A `-s` override shows up in the outputs.
    - 13 new loader tests; 45 tests in total.
- 2026-10-23, **Phase 3 done** (details in the report, Phase 3 section)
  - **Input checks:** `validate_inputs.py` holds the profile and climatology contracts
    (`docs/input_contract.md`).
    - `NCEI.py` checks the climatology, then all input files, before processing; failing files are
      skipped with reasons.
    - Each worker returns written / no data left / failed, with the step and error.
    - `<run>/documentation/run_summary.txt`; exit code 3 when files failed (`run_ncei.sh` still runs
      steps 11–13 on the good files).
  - **Climatologies:**
    - Any file meeting the contract can be used; `setup_climatology.py` checks, builds (`--build`) or
      pre-bakes (`--prebake`).
    - A missing climatology stops the run with instructions.
    - Pre-baked files are used only if their `source_file` matches the current climatology. Before,
      a stale pre-baked file could silently be used.
  - **Cleanup:** pdb removed; tracebacks logged.
  - **Found in the real inputs:** 18,374 CTD profiles (1.6%; 1992–2002, up to 7.3% in 2002) with the day missing (dates
    `YYYYMM00`), zero-weighted by step 06 as before. Flagged for the team; not changed.
  - The user's decision (2026-10-23): note the mid-month option (the 15th), but DON'T implement it
    yet.
  - **Tests:** 76 in total.
    - 27 contract tests;
    - 3 climatology-guard tests;
    - a bad-file end-to-end test;
    - the regression test still matches the Phase 1 reference (`run_summary.txt` added).
  - **Deferred to Phase 4:** print → logging. The log format is parsed by
    `diagnostics/lib/costs_scaling.py`, so keep its lines.
- 2026-10-23, follow-up to Phase 3 (from the user)
  - A comment box above `paths.climatology_file` in `config/defaults.yaml`, `local.example.yaml` and
    `local.yaml`: use your own file, or `setup_climatology.py --build`.
  - A new **Quick start** near the top of `README.md`, with "you don't need a climatology already" as
    step 2.
  - The Sphinx front page must repeat this (see Phase 5).
- 2026-10-23, **Phase 4 done** (details in the report, Phase 4 section)
  - **Package:** the code is now `ecco_insitu/` (`pipeline`, `steps/`, `packing`, `config` + `defaults.yaml`,
    `validate`, `climatology/`, `diagnostics/`, `cli`, `log`), moved with `git mv`; imports
    rewritten; no `sys.path` hacks.
  - **Install and command:** `pyproject.toml` (deps + extras); `pip install -e .`; the `ecco-insitu` command
    with 16 subcommands. 17 old script names are kept as wrappers. `run_ncei.sh` runs
    `python -m ecco_insitu...`.
  - **Logging:** controller messages go through logging (`ecco_insitu/log.py`). The steps' per-file output
    is still printed verbatim, because `costs_scaling.parse_step_counts` parses it; checked
    identical before and after.
  - **Tests:** 111, regression bit-identical to the Phase 1 reference; `ecco-insitu run` from `/tmp`
    is identical too.
  - Installed editable in `ecco_env1`.
- 2026-10-23, **Phase 4 follow-up** (the user is the sole user and asked for the modern design, no shims)
  - **Removed:** the 17 wrapper scripts and `run_ncei.sh`.
  - **New:** `src/ecco_insitu/run.py` (`ecco-insitu run`), with OS-level tee to the
    `z_logs/ncei_log_<ts>.log` (`run.log_dir`) and the same exit codes; flags renamed (`--no-scaling`,
    `--plots`, `--compact`, `--no-pack`).
  - **Moved:** the `src/` layout; `diagnostics/README.md` → `docs/diagnostics.md`.
  - **Fixed:** a hang when the run finished (the resource tracker held the tee pipe), now an end marker.
  - **Tests:** 95 pass, still bit-identical to the Phase 1 reference.
  - **Stale:** `patches/step06_code12...` no longer applies (noted in `HANDOFF.md`).
- 2026-10-23, **Phase 5 done** (details in the report, Phase 5 section)
  - **Docstrings:** every module, public function and class (4 parallel workers). Verified docstring-only
    with `tests/tools/same_code.py docstrings-base` (39 files, 0 code changes).
  - **Sphinx:** in `docs/sphinx/`, rtd theme, napoleon, myst.
    - The front page has the "you don't need a climatology" box.
    - The user guide, plus an API reference of 43 module pages.
    - Builds with 0 warnings.
  - **New:** `CHANGELOG.md`; `tests/unit/test_docs_build.py` (docstring coverage + strict build). 97 tests.
  - **Issues found:** listed in report section 5.4, none fixed. Safe fixes: `np.NaN` (csv_to_nc,
    LLC270 loaders), the `logical_and` bug, stale texts, unused imports. Team decisions: csv_to_nc
    T label and years, step06 T ≥ −2 °C and (0,0), step09 HHMMSS distance.
- 2026-10-23, **known issues published** (the user's request)
  - `docs/known_issues.md` (KI-1 to KI-17) is also a docs page ("Known issues").
- 2026-10-23, **Phase 6 done** (details in the report, Phase 6 section)
  - **CI:** `.github/workflows/tests.yml` runs on push/PR, py3.11 and 3.12; installs with
    `.[diagnostics,test,docs]`, runs pytest (strict docs build included) and `ecco-insitu --help`.
    - Data tests skip there (5); checked by simulation: 92 passed, 5 skipped, also in a fresh venv.
  - **Docs publishing:** `.github/workflows/docs.yml` publishes to GitHub Pages from `main`. It needs the
    repo setting Pages → Source: GitHub Actions, which only the owner can set.
  - **PR checklist:** `.github/pull_request_template.md` (science / changelog / docstrings / docs /
    known issues).
  - **Version stamp:** `ecco_insitu_version` global attribute in every output NetCDF, plus the
    `run_config.yaml` header.
    - Reference regenerated deliberately: data identical; only `ecco_insitu_version` (new) and
      `created_by` changed.
  - **Tests:** find the data via the settings (no `/Users` paths in tests).
- **Project complete through Phase 6.** Suggested next steps:
  - merge to `main` (the user's call), enable Pages, tag `v1.0.0`;
  - the safe known-issue fixes (KI-1, 2, 10, 11, 13–17);
  - team decisions on KI-4–9;
  - optionally, trimmed data so the regression tests run in CI.
