# Multi-source Profile Globe App

Interactive Dash/Plotly application for exploring ECCO in-situ profile data
across instrument sources, geodesic bins, and time.

## What the app shows

- A rotatable orthographic globe with geodesic bins coloured by sampling
  probability (fraction of months in the dataset record that received at least N
  profiles in that bin).
- Two resolution modes: **fine** (10242-bin geodesic grid, `prof_bin_id_a`) and
  **coarse** (2562-bin grid, `prof_bin_id_b`).
- Clicking a bin loads depth–time heatmaps of monthly-mean T or S (observed,
  climatology, or anomaly) for that bin, plus an optional source-profile table.
- Data is partitioned by instrument source (CTD\_WOD, GLD\_WOD, MEOP, MRB\_WOD,
  PFL, PFL\_BGC, Samoa, XBT\_WOD). The "All (combined)" view pools sources with
  profile-weighted averaging.

---

## Files

```
profile_analysis/
  app_allsource_globe.py          — the app itself
  build_allsource_bin_timeseries.py  — preprocessing script that builds the pickle
  README_globe_app.md             — this file
```

---

## Dependencies

### Python packages

```
dash
plotly
numpy
scipy          (SphericalVoronoi for geodesic cell polygons)
pandas
xarray         (build script only)
```

Install with:

```bash
pip install dash plotly numpy scipy pandas xarray
```

### Data files

Two external data files are required. **Neither is included in this repository.**

| File | Purpose | Approx. size |
|------|---------|-------------|
| `allsource_bin_timeseries.pkl` | Pre-aggregated profile data (all sources, both grids, monthly means) | ~4 GB |
| `02562_bin_locations.csv` | Lon/lat centres for all 2562 coarse-grid geodesic bins | <1 MB |

The geodesic CSV files (`10242_bin_locations.csv`, `02562_bin_locations.csv`)
are part of the ECCO in-situ support data package (available on the project
Google Drive).

The pickle is produced by running `build_allsource_bin_timeseries.py` (see
below). It is not a raw data file; it is a derived product from the NCEI
post-processing chain output.

---

## Configuring paths

The app has two hardcoded paths near the top of `app_allsource_globe.py`:

```python
PICKLE_FILE       = Path('/path/to/allsource_bin_timeseries.pkl')
GEODESIC_2562_CSV = Path('/path/to/02562_bin_locations.csv')
```

Edit these before running. The rest of the app reads everything from the pickle
and does not require any other files at runtime.

---

## Running the app

```bash
cd profile_analysis
python app_allsource_globe.py
```

Then open `http://127.0.0.1:8050` in a browser.

Startup takes 30–60 seconds: the pickle (~4 GB) must load and two
SphericalVoronoi tessellations (10242 + 2562 cells) must be computed before
the server accepts connections.

---

## Building the pickle from scratch

If you have your own NCEI-processed profile files and want to build the pickle:

```bash
python build_allsource_bin_timeseries.py \
    --source_root /path/to/profile_files_NCEI_processed \
    --output      /path/to/allsource_bin_timeseries.pkl
```

The build script also requires two geodesic CSV files. Their paths are
currently hardcoded in the script (`_DEFAULT_GEODESIC_FILE`,
`_DEFAULT_GEODESIC_FILE_B`) and must be edited before running.

### Input: NCEI-processed profile files

The source directory must contain one subdirectory per instrument source
(e.g. `CTD_WOD/`, `PFL/`), each holding NetCDF files produced by the NCEI
post-processing chain (`NCEI.py`). Each file must contain at minimum:

- `prof_YYYYMMDD`, `prof_lon`, `prof_lat`
- `prof_T`, `prof_Tclim` (and `prof_S`, `prof_Sclim` where available)
- `prof_bin_id_a` — 10242-bin geodesic bin index (assigned by NCEI chain step 2)
- `prof_bin_id_b` — 2562-bin geodesic bin index (assigned by NCEI chain step 2)

Producing these files from raw WOD data requires running the full NCEI
processing chain (`NCEI.py` and its 10 steps), which depends on additional
support files (LLC90 grid, WOA23 climatology, sigma files). See the top-level
`README.md` for details on that pipeline.

### Output

The build script writes a single pickle containing:

- Per-bin, per-source, per-month mean profiles of T, S, T\_clim, S\_clim
- Long-run T/S means and sampling hit counts per bin
- Both fine (10242) and coarse (2562) grid aggregations
- Dataset date range and geodesic bin centre coordinates

Build time is approximately 15–20 minutes for ~420 files across 8 sources.

---

## Controls

| Control | Effect |
|---------|--------|
| Drag globe | Rotate |
| Scroll wheel | Zoom toward cursor |
| Zoom slider | Zoom (alternative) |
| Grid resolution | Switch between 10242-bin (fine) and 2562-bin (coarse) geodesic grids |
| Source | Filter to a single instrument source or show all combined |
| Globe: probability threshold | Colour bins by P(≥1), P(≥2), or P(≥3) profiles/month |
| Probability basis | Denominator: full dataset span, or since each bin's first observation |
| Field | Temperature or salinity |
| Anomaly reference | Subtract climatology or bin long-run time-mean |
| Click a bin | Load depth–time heatmaps and (optionally) source-profile table |
| Show source table | Replace anomaly heatmap with a table of individual profiles |

---

## Notes on data coverage

- Bins with no data for the selected source are not shown on the globe.
- Switching grid resolution clears the current bin selection.
- The coarse grid (2562 bins) will show higher probabilities than the fine grid
  (10242 bins) because each coarse bin aggregates the profiles of roughly four
  fine bins.
- Model-equivalent profiles (`prof_Testim`, `prof_Sestim`) are not yet
  populated in the upstream NCEI chain and are not displayed.
