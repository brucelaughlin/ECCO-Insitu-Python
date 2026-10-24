Automatic checks and publishing
===============================

The repository has two automatic jobs, defined in ``.github/workflows/``. They run on GitHub's
servers, not on your computer. GitHub reads these files when commits are pushed, and runs them
itself. Nobody has to start them.

How it works
------------

`GitHub Actions <https://docs.github.com/actions>`_ is a service built into GitHub. A *workflow*
file says **when** to run (for example, "on every push") and **what** to run (a list of shell
commands). When a matching event happens, GitHub:

1. starts a fresh, empty virtual machine (Ubuntu Linux) in its data center;
2. copies the repository at that commit into it;
3. runs the workflow's commands;
4. records the result (passed / failed) and the full log, then deletes the machine.

This is free for public repositories. The results are on the repository's **Actions** tab, and as
a green tick or red cross next to each commit and pull request.

Who runs it? Nobody, directly
-----------------------------

Nothing on your computer, and no person, "controls" these runs:

- **Starting a run:** pushing a commit is the only action. When the push includes a workflow file, or
  changes the code, GitHub runs the workflows by itself on its own machines. Nobody (including an
  assistant editing the code) logs into those machines or starts anything there.
- **Following a run:** watch it on the repository's **Actions** tab, which shows each run, each step
  and the full log. The same information is in GitHub's public API for a public repository.
- **When it runs:** every push to any branch triggers the ``tests`` workflow, whoever pushes. The
  ``docs`` workflow runs only when ``main`` changes (or by hand), and publishes only once GitHub
  Pages is enabled.

The first run (2026-10-23) followed the push that added the workflow files. Before that push, the
same checks were simulated locally, so the first run wouldn't fail: a clean copy of the repository
with no data, and a brand-new Python environment. It passed on Python 3.11 and 3.12.

What it runs, and what it does not
----------------------------------

The automatic run is **more than installing the package, but it does not run the NCEI chain on
profiles**:

- **Installs the package from scratch.** This proves ``pyproject.toml`` lists everything the package
  needs.
- **Runs the test suite** (92 tests there), which exercises the code:

  - the input checker, on deliberately broken files (missing variables, reversed depths,
    temperatures in kelvin, ...);
  - the settings loader (typos, wrong types, merge order);
  - small unit tests on individual pieces: date/time rollover, step 05's potential-temperature
    decision, the UNESCO check value of the potential-temperature formula, step 11's time windows,
    and every ``ecco-insitu`` command starting;
  - the docstring check and the strict documentation build.

- **Starts** ``ecco-insitu --help``.

The two tests that run **the whole chain** on the 8 small test files are *skipped* there. They need
the grid, climatology and σ files (about 1.5 GB), which are not in the repository. So:

- **GitHub checks**, on every push, that the package installs, its pieces behave, and the docs
  build.
- **Only a machine with the data checks** that the chain still produces bit-identical scientific
  output: ``python -m pytest`` with ``config/local.yaml`` pointing at the data.

**Possible later step:** commit trimmed copies of the data, cut down to what the 8 test files touch
(tens of MB rather than 1.5 GB), so the full-chain tests run on GitHub too.

The ``tests`` workflow (``.github/workflows/tests.yml``)
---------------------------------------------------------

**When:** every push to any branch, and every pull request.

**What:** on Python 3.11 and 3.12 separately:

- install the package as a new user would: ``pip install -e ".[diagnostics,test,docs]"``;
- run ``python -m pytest``, including:

  - the unit tests;
  - the check that every public function has a docstring;
  - the documentation build in strict mode (any warning fails);

- run ``ecco-insitu --help``.

**What it can't check:** the tests that need the chain's external data (grid, climatology, σ files,
about 1.5 GB, not in the repository) are *skipped* there, with a reason in the log:

- the two full-chain regression tests;
- the checks of the real climatology.

They run on any machine whose ``config/local.yaml`` points at the data. So **run**
``python -m pytest`` **locally before merging** a change that could affect results; the
pull-request checklist asks for this.

The ``docs`` workflow (``.github/workflows/docs.yml``)
-------------------------------------------------------

**When:** every push to the ``main`` branch, or by hand (Actions tab → *docs* → *Run workflow*).

**What:** builds this documentation (strict mode) and publishes it to GitHub Pages:
https://brucelaughlin.github.io/ECCO-Insitu-Python/

**One-time setup** (repository owner): *Settings → Pages → Source: GitHub Actions*.

The pull-request checklist
--------------------------

``.github/pull_request_template.md`` fills in the description of every new pull request with a
checklist:

- the tests pass locally, with the data;
- whether the change alters the processed profiles, and if so, that the regression reference was
  regenerated on purpose and the reason recorded in ``CHANGELOG.md``;
- a changelog entry;
- docstrings;
- docs;
- known issues.

How the pieces fit
------------------

==========================================  ===================================================
Kind of mistake                              What catches it
==========================================  ===================================================
a change alters processed values             regression test (locally, with the data)
a public function has no docstring           test, locally and on GitHub
a docstring breaks the documentation         strict docs build, locally and on GitHub
a missing dependency                         the clean install on GitHub
a forgotten changelog / docs update          the pull-request checklist
==========================================  ===================================================

Running the same checks yourself
--------------------------------

.. code-block:: bash

   python -m pytest                                                # everything (~90 s, with data)
   python -m sphinx -b html -W docs/sphinx docs/sphinx/_build/html  # the docs, as on GitHub
