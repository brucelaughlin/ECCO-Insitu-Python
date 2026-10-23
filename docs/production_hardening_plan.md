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
- **Next: Phase 2 (configuration).**
  - Move every path and option out of `NCEI.py` and `run_ncei.sh` into one YAML file, with defaults
    equal to today's values, so the regression test still passes unchanged.
  - Save the config used into each run's output.
  - Make the data paths environment-overridable, so CI and other users can point them elsewhere.
