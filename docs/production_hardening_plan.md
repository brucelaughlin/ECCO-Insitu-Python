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
  - **Found in the real inputs:** about 18,000 CTD profiles (1992–2002) with month-only dates
    `YYYYMM00`, zero-weighted by step 06 as before. Flagged for the team; not changed.
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
- **Next: Phase 4 (package).**
  - Move into `ecco_insitu/` (`steps/`, `config`, `validate`, `io`, `climatology`, ...), with a
    `pyproject.toml` and an `ecco-insitu` command line (run / check / setup-climatology).
  - Keep `run_ncei.sh` and the old script names working as thin wrappers for now.
  - Switch to `logging` while keeping the log lines the diagnostics parse.
  - The regression test must pass unchanged.
