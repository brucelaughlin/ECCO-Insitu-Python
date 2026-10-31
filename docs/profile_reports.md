# Profile reports: how they are made and published

This document describes the profile analysis reports (the interactive web pages on GitHub Pages) as of
2026-10-30: what makes them, what is in them, the choices made, and how publishing through GitHub works.
It is not part of the documentation website.

Code: `src/ecco_insitu/diagnostics/profile_analyses/` (`compute.py`, `report.py`, `viewer.py`,
`viewer.html`, `publish.py`, `compare.py`).

## 1. From a chain run to a published report

```bash
# 1. model equivalents (approximate, from v4r6 monthly means), once per run
ecco-insitu monthly-estim --theta <THETA_mon_mean dir> --salt <SALT_mon_mean dir> \
    --data <run>/ecco_arrays --out <equiv dir>/netCDF

# 2. the analyses: reads every profile once per pass (two passes), ~40 min for a full run
ecco-insitu profile-analyze --run <run>/yearly --model-equivalents <equiv dir>/netCDF --out <analysis dir>

# 2b. figures only, from the saved results (minutes): after any change to labels, layout or colors
ecco-insitu profile-analyze --replot --out <analysis dir>

# 3. the web pages
ecco-insitu profile-report --site --inputs <analysis dir> --out docs/reports/<report name>

# 4. publish (section 4)
ecco-insitu publish-reports docs/reports/<report name> [docs/reports/<other report> ...]
#    then on GitHub: Actions -> docs -> Run workflow (branch main), or push anything to main
```

A full re-analysis (step 2) is needed only when a calculation changes. Everything else (titles, legends,
colors, layout, which split the core figures use) is a replot (2b) plus steps 3–4.

Where things are on this machine (2026-10-30):

| What | Where |
|---|---|
| analysis results (NetCDF + figures) | `~/ecco/yip/insitu_reports/profile_analyses_20261028_142247/` (run 142247), `.../profile_analyses_20261027/` (run 011321) |
| report sites (made by step 3) | `ECCO-Insitu-Python/docs/reports/profile_report_20261027_142247/`, `.../profile_report_20261027/` |
| model equivalents | `~/ecco/yip/ecco_v4r6_model_data/model_equivalents_20261028_142247/`, `.../model_equivalents_20261027/` |

## 2. What `profile-analyze` computes

Per obs class (CTD, glider, ITP, MEOP, mooring, Argo, XBT) and cost type (`vs_clim`: WOA23 1991–2020;
`vs_model`: v4r6 monthly means sampled at the profiles; `vs_run_b`: a second set of monthly means), for
T and S:

- **Weight versions:** `final` (γ·sf/σ², as ECCO), `nosf` (γ/σ²), `true` (1/σ², "cost per datum").
- **Unsplit** cost: w·(V − V_ref)².
- **Split** cost, as ECCO v4r6's `cost_profiles.F`, in two terms, with V̄ and V̄_ref the time means of the
  points in each geodesic bin and depth level (all years pooled):
  - **anomaly term**, per point: w·((V − V̄) − (V_ref − V̄_ref))²;
  - **mean term**, per bin and level: w̄·(V̄ − V̄_ref)², w̄ = 1/max(s/n, s/9), s = mean(1/w), n = points.
- The split is computed twice: with each class's **own** bin means, and with bin means **pooled** over all
  classes (what `cost_profiles.F` does). Which one the core figures show is a setting (section 3).
- Everything is saved in one NetCDF per class and cost type (`analysis_<class>_<run>_<cost type>.nc`),
  including the bin means, the mean term per (bin, level), cost per bin and depth band, and, for the map
  explorer, the cost per datum and the signed misfit summed per non-empty (year, bin, level), stored
  sparsely (`viewer__*` variables).

**Checked:** the pooled anomaly cost and mean term equal `ecco-insitu model-cost` (which reproduces the
team's v4r6 cost bars) on run 20261027 vs_model, to the reported precision; a unit test checks the same
on synthetic data. Not yet checked against a model run's own cost output.

## 3. The report and the choices made

The report page (`index.html`): tabs for obs class, cost type, and analysis (01–10 the team's core list;
11 and 12 supplemental), with extra rows where they apply (**Cost**: unsplit / anomaly / mean term for 06
and 07; **Depth**: standard depth bands for 05 and 06). Figures share axis and color ranges across cost
types, and have the same layout in every class. Click a figure to see it full width. The side panel has
the description and Terminology; **Notes & errata** (top right) has interpretation, data notes, open
decisions and corrections; **Map explorer** (top right) is the interactive map page (below).

Choices made, all reversible:

| Choice | What we did | To change |
|---|---|---|
| Split in the core figures (2026-10-29) | each class's **own** bin means, so both terms are per class (the team's per-class rule); the pooled `cost_profiles.F` version is analysis 12 and, for the mean term, in Metadata | `ecco-insitu profile-analyze --replot --split pooled --out <dir>` (or `CORE_SPLIT` in `compute.py`); no re-analysis |
| Depth bands for the map images (2026-10-29) | 0–50, 50–100, 100–250, 250–500, 500–1000, 1000–2000, 2000–4000 m and deeper, plus all depths | `MAP_DEPTH_BANDS` in `compute.py` (re-analysis) |
| Map explorer time resolution (2026-10-31) | 1-year blocks, 1992–2026 (was 5-year until 2026-10-31); stored sparsely | `VIEWER_YEAR0`, `VIEWER_BLOCK_YEARS`, `VIEWER_N_BLOCKS` in `compute.py` (re-analysis) |
| Mean-bias diagnostic (2026-10-31) | signed misfit per bin, any years; shown separately from the costs | `_add_viewer` in `compute.py`, `viewer.html` |
| Flagged spikes (2026-10-30) | a level whose mean cost per datum (true; unsplit or anomaly; ≥ 30 points) is > 0.5 and ≥ 1.4× the median of the levels 300–1000 m above and below | `_spikes` in `viewer.py` (rebuild the site only) |
| Analysis 1 (2026-10-28) | every weighted observation, the same in every cost-type tab | see Notes & errata in the report for the old behavior and how to revert |

## 4. The map explorer

`viewer.html` draws the maps in the browser instead of loading images, so any combination is available:
obs class, cost type, variable; cost (unsplit, split anomaly, split mean term, number of points); depth
(all, any single level, 100 m bands, the 97 levels in 15 or 30 groups, or the automatically flagged spike
levels); time (all years, or any range of years with the two-handled slider, or one year at a time; not
for the mean term, which is defined over the whole record). Hover a bin for its value; the color bar has
ticks; a log color scale and a fixed maximum are available. A **Guide** on the page explains every option.

**Mean-bias diagnostic** (its own row, separate from the costs): the signed misfit, observation minus
reference, averaged per bin over the selected depths and years, in °C or psu, on a red–blue scale centered
on zero. It is not a cost term. ECCO's mean term is defined over the whole assimilation period (the window
the estimate is fit over), so a per-period mean *term* would not be ECCO's; but the signed offset over
chosen years shows whether it drifts (model drift, instrument or region changes, ocean change against a
fixed climatology). Checked: over all years it reproduces the per-bin bias saved by `profile-analyze`.

How it works: `profile-report --site` writes `viewer/`:

- `grid.json` and `grid_<view>.bin`: for each view (global equal area from Australia's western tip,
  north and south polar stereographic), the geodesic bin under each pixel, and the coastline;
- `<class>__<cost type>.bin`: for every non-empty (bin, level, year) the mean cost per datum unsplit and
  anomaly, the point count and the mean signed misfit; and the mean term per (bin, level). Stored as gzip of 16-bit
  values (3 significant digits; values below ~6×10⁻⁵ lose precision but are zero on every color scale).
- `index.json`: the files, depth levels, blocks, level groups and flagged spikes.

The browser downloads one class's file when it is first selected (sizes in section 4 of the HANDOFF /
the build output: a few hundred kB for small classes, tens of MB for Argo). **Checked:** decoding the
files reproduces analysis 02's mean cost per datum at the levels tested.

The explorer needs a web server (GitHub Pages, or `python -m http.server` in `docs/reports/` for a local
look); opening `viewer.html` straight from the disk doesn't work, because browsers block file downloads
from `file://` pages.

## 5. Publishing through GitHub

### The pieces

- **The repository** `brucelaughlin/ECCO-Insitu-Python` (GitHub). Branches are separate lines of commits.
- **`main`**: the code and the documentation sources. Everything is developed here.
- **`gh-pages`**: only the published report sites. It has **no history**: every publish replaces it with
  one new commit (a force push), so it always holds exactly one copy of the reports.
- **GitHub Pages**: GitHub's web hosting. It serves whatever a deploy uploads, at
  `https://brucelaughlin.github.io/ECCO-Insitu-Python/`.
- **GitHub Actions**: GitHub's automation. The workflow `.github/workflows/docs.yml` ("docs") is what
  deploys to Pages. The repository's settings (Settings → Environments → github-pages) allow deploys
  from `main` only, so the workflow always runs on `main`.

### What the docs workflow does

1. Checks out `main` and builds the documentation website (Sphinx, from `docs/sphinx/`).
2. If the `gh-pages` branch exists, copies its report folders into the built site.
3. Uploads the result to GitHub Pages.

It runs on every push to `main`, and when started by hand (Actions → docs → Run workflow).

So the documentation always comes from `main` and the reports always come from `gh-pages`, and any deploy
publishes both.

### Why the reports are not on `main` any more

Git keeps every version of every file on a branch forever. The report sites are large (~170 MB each,
mostly figures) and change completely at each rebuild, so committing them to `main` (as until 2026-10-29)
added their whole size to the history every time: the repository grew to ~210 MB in three days. On
`gh-pages`, which is replaced at each publish, old versions are discarded and nothing accumulates.

`docs/reports/` is now in `.gitignore`: the folder still exists on this machine (it is where
`profile-report --site` writes), but git on `main` ignores it.

The old report images are still in `main`'s history (~210 MB). Removing them would mean rewriting the
history of `main` (every commit changes, and every clone has to be replaced); not done.

### Day to day

| To | Do |
|---|---|
| update a report | `profile-report --site ... --out docs/reports/<name>`, then `ecco-insitu publish-reports docs/reports/<name> [other reports]`, then Actions → docs → Run workflow (or push to `main`) |
| change code or documentation | commit and push `main` as usual; the deploy keeps the published reports |
| add a report (new URL) | as above with a new `<name>`; list every report to keep, or use `--keep-others` |
| remove a report | publish the others without it (or delete it locally and publish what is left) |
| check before publishing | `ecco-insitu publish-reports ... --dry-run`, or look locally: `cd docs/reports && python -m http.server`, then open `http://localhost:8000/<name>/` |

A report's URL is `https://brucelaughlin.github.io/ECCO-Insitu-Python/<name>/`. Figure links carry a
version tag from the image's content (`...png?v=<hash>`), so a changed figure is never shown from a
browser's cache.

`publish-reports` starts the deploy itself if the GitHub command-line tool `gh` is installed and logged in
(`brew install gh`, `gh auth login`); otherwise it prints the one manual step.

### If something goes wrong

- **A report URL gives 404:** the deploy hasn't run since the publish (run the docs workflow), or the
  report wasn't in the last publish (publish it again, with the others).
- **The deploy failed** (red X in Actions): open the run to see the step. A failure in "Build" is the
  documentation (Sphinx treats warnings as errors); the reports are only copied in afterwards.
- **Old figures show:** reload; the version tags make this rare.
