# Profile map app

`app_profile_map.py` is an interactive Dash/Plotly app for exploring the NCEI-processed
in-situ profiles by source, geodesic bin and time, on one flat map. It is **standalone**.
It reads a data store (one folder: zarr + JSON + SQLite) that
`build_allsource_bin_timeseries_with_anomalies.py` builds from post-NCEI profile files.
Anyone with those files, or with a ready store, can run it without the NCEI chain: see
**[What you need](#what-you-need)**.

It replaces the three earlier globe apps (`app_allsource_globe.py`,
`app_allsource_anomaly_globe.py`, `app_profile_scatter.py`). They are kept at git tag
`globe-apps-2026-10-19`; restore them with
`git checkout globe-apps-2026-10-19 -- profile_analysis/`.

## What it shows

- **Projection:**
  - **World:** natural earth.
  - **Arctic** and **Southern Ocean:** polar stereographic, centered on the pole. They
    open showing poleward of 50° (`map_common.POLAR_LAT_LIMIT`).
  
  Plotly draws these projections in a square area, so the polar maps sit in the middle of
  the wide map panel. Changing projection resets the zoom.
- **Show:**
  - **Bins:** geodesic bins colored by one of:
    - visit probability P(≥N), described below;
    - the long-run T or S anomaly (value − climatology) **mean** or **std** at a chosen
      depth, from the profiles of the **selected source** (or all sources pooled). A bin
      where that source has no data at the chosen depth is left unfilled.
    
    In the **anomaly-mean** modes two more things are drawn (described below): a
    **std outline** inside each bin, and **sub-polygons** when zoomed in.
  - **Profiles:** individual profiles as dots colored by source, within a year window,
    over faint bin outlines. About 20,000 dots are drawn, thinned evenly over **the area
    in view**, and re-chosen as you pan and zoom, so a zoomed-in area shows its own
    profiles densely. The header says how many of the profiles in view are shown.
- **Click a bin or a dot:** the panels below show that bin's depth–time T / S
  (observed or climatology) and its anomaly (− climatology or − bin time-mean). A
  table of the bin's individual profiles is also available.
- **Two grids:** fine (10242 bins, `prof_bin_id_a`) and coarse (2562 bins, `prof_bin_id_b`).
- **Sources:** one per source folder of the run's `yearly/` directory (today CTD_WOD,
  GLD_WOD, ITP, MEOP, MRB_WOD, PFL, XBT_WOD). Files in sub-folders (ITP/L2, ITP/L3)
  belong to their top folder's source. "All (combined)" pools the sources, weighting
  each by its number of profiles.

### Visit probability

P(≥N) is the fraction of months in which a bin had at least N profiles:

> P(≥N) = (months in which the bin had at least N profiles) ÷ (months considered)

- **N** is 1, 2 or 3 ("a month counts if the bin has at least N profiles"). P(≥1) is the
  chance that a given month has any profile in the bin.
- **Months considered** ("out of"):
  - **all months of the record:** January 1992 to August 2026 (416 months) for every bin;
  - **the months since the bin's first profile:** this doesn't penalize bins, or
    sources, that started being sampled later (e.g. Argo from the late 1990s).
- **Source:** with a source selected, only its profiles count, and "since the first
  profile" starts at that source's first profile in the bin.
- **Depth:** a profile counts for its month whatever depths it covers; there is no depth
  choice in this mode.
- **Grid:** coarse bins cover about four fine bins, so they show higher probabilities.
- **Example:** P(≥1) = 0.8 means at least one profile in 80% of those months.
- **In the app:** a short explanation under the probability controls restates this for the
  current choices. The probability also appears in the hover and bin info in the anomaly
  modes.

### Inset outlines (anomaly modes)

Each bin has an outline drawn 8% inside its edge (`EDGE_INSET`), so neighbouring bins never
share an outline. The note at the bottom left of the map says what they show:
- **Width (zoomed out):** proportional to the number of the selected source's profiles
  with data at the chosen depth, from 0.2 px to 3 px (reached by the busiest 2% of bins),
  so densely sampled areas stand out. Widths grow by √2 per zoom doubling.
- **Width (sub-polygons showing):** a uniform 1.5 px, since only the color matters then.
- **Color:** in the anomaly-mean modes, the bin's anomaly **std** at that depth (cividis,
  0 to the 98th percentile; its color bar is on the left of the map, the fill's on the right). In the
  std modes the fill already is the std, so the outlines are light gray.
- **No data at this depth:** bins where the selected source has profiles, but none at
  the chosen depth, are left unfilled with a thin gray outline, so the source's footprint
  stays visible. The note says how many bins have data at that depth (e.g. ITP starts a
  few meters down, so at 2 m it has data in only 1 of its 110 bins).
- **Turn them off** with "Edges: anomaly std". The outlines stay on top of the
  sub-polygons.

### Sub-polygons (anomaly-mean modes, zoomed in)

**When they appear:** once the view is less than 15° of latitude tall
(`SUBPOLY_VIEW_DEG`), the same for every source and projection: about zoom 12 on the
world map and 5.6 on a polar map. At that size even the densest area has at most ~150
bins in view. `SUBPOLY_MAX_BINS` = 150 stays as a safety limit; if it is hit, the
header asks you to zoom in a little more.

**What they show:** each bin is split into cells showing the spread of its single
profiles at the chosen depth, as in the original `geodesic_binning/` tool:
- **n ≤ 100 profiles:** one cell per profile.
- **More profiles:** 100 cells, each one quantile ((i + 0.5)/100) of the bin's anomalies.
- **Source:** like the bin colors and outlines, the cells use only the selected
  source's profiles (or all sources pooled).
- **Order:** values run low to high from the south of the bin to the north (near the
  poles: from the edge toward the pole).
- **Hover:** a cell shows its bin and depth, which profile it is (e.g. "Profile 3 of 52", counted from the
  lowest anomaly), or which percentile of
  the bin's profiles it stands for ("Percentile 6.5 of 1,234 profiles"), and
  its value as T (or S) − climatology.

**Color scale:** single profiles vary far more than bin means, so the cells get their
own symmetric range from the cells in view. Their color bar ("single profiles") replaces
the bins' bar while they show.

**Geometry:** Voronoi cells with Lloyd relaxation (`subpolygons.py`), computed in a
tangent-plane projection around each bin, so polar and dateline bins work too. They are
computed when first needed (~4 s for 100+ bins) and cached (~0.1 s afterwards).

---

## What you need

The app is **standalone**. It needs three Python files from this repo and a **data store**
(one folder). The store is built from **post-NCEI profile files**; you don't need to run
the NCEI chain yourself if someone gives you its output. Start from whichever you have:

| You have | Then |
|---|---|
| **A data store** (a folder holding `<name>.zarr`, `<name>_metadata.json`, `<name>_profiles.db`), e.g. copied from a colleague | Run the app: step 3 |
| **Post-NCEI profile files**, from a chain run (yours or someone else's) | Build a store (step 2), then run the app |
| **Only raw profile files** (UCSD `Interp_Profiles`) | Run the NCEI chain first (top-level `README.md`), then steps 2 and 3 |

### 1. Files from this repo

- **To run the app:** `profile_analysis/app_profile_map.py`, `map_common.py` and
  `subpolygons.py`. Nothing else from the repo is needed: not the chain, `tools.py` or the
  model grid.
- **To build a store:** also `profile_analysis/build_allsource_bin_timeseries_with_anomalies.py`.
- **Python packages:**
  - app: `dash`, `plotly`, `numpy`, `scipy`, `zarr` (3.x), `shapely` (2.x), `matplotlib`
    and `cartopy`;
  - build: also `xarray`, `netCDF4` and `pandas`.
  
  Tested with zarr 3.1.5, dash 4.4.1, plotly 7.0.0, shapely 2.1.2, cartopy 0.25.0 and
  scipy 1.17.1.

### 2. Build a store from post-NCEI files

**Inputs:**
- **The post-NCEI profile files** (`*.nc`), in either layout:
  - a chain run folder `profile_files_NCEI_processed_<ts>/`, whose `yearly/` is used
    (`ecco_arrays/` and `documentation/` are not needed);
  - or any folder with **one sub-folder per source**, `<folder>/<SOURCE>/*.nc`. Files may sit
    deeper inside a source (e.g. `ITP/L2`, `ITP/L3`), and the sub-folder name is what the
    app shows as the source.
- **What each file must contain:**
  - `prof_YYYYMMDD`, `prof_lon`, `prof_lat`, `depth` (the 97-level grid, 2–6000 m);
  - `prof_T`, `prof_Tclim`, `prof_S`, `prof_Sclim` (a file may lack T or S);
  - `prof_bin_id_a` / `prof_bin_id_b`, the geodesic bins.
  
  The climatology and the bin numbers are added by the NCEI chain (steps 02–03). That is
  why the store needs post-NCEI files, not raw ones.
- **The two geodesic bin-center files,** `10242_bin_locations.csv` and
  `02562_bin_locations.csv`, from the ECCO in-situ support data (the `geodesic` folder on
  the project Google Drive). Pass their folder with `--geodesic_dir`.

**Command:**
```bash
python profile_analysis/build_allsource_bin_timeseries_with_anomalies.py \
    --source_root <run folder, or folder of source folders> \
    --geodesic_dir <folder with the two bin-center files> \
    --output <store folder to write>
```

- **Output:** a store folder holding `<name>.zarr`, `<name>_metadata.json` and
  `<name>_profiles.db`.
  - Default name: `allsource_bin_timeseries_with_anomalies_<ts>`, where `<ts>` is the
    run's timestamp.
  - Default location: this project's store folder if it exists, else the current folder.
- **Cost (for one full run, ~6 million profiles):** about 25–30 min and a 7.4 GB store.
  The build holds every profile in memory; it ran comfortably on a 32 GB machine.
- **Rebuild after each chain run you want to look at;** a store always shows the run it
  was built from.
- **Adding to an existing store:** `--append_anomalies` adds or refreshes the anomaly
  arrays in place, re-reading the profile files, without rebuilding the rest (used when
  new arrays were added to the format).
- **If a build hangs** after "Aggregating N bins...", rerun it with
  `OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1`.
- The script stops with a clear message if `--source_root` or a bin-center file is missing.

### 3. Run the app

```bash
python profile_analysis/app_profile_map.py --data <store folder>
```

Then open `http://127.0.0.1:8050` in a browser.
- **Default store:** without `--data`, it loads the newest
  `allsource_bin_timeseries_with_anomalies_<ts>` store in this project's store folder
  (`/Users/brucel/ecco/yip/profile_data/z_profile_file_analysis/`). It says so if there is
  none.
- **First run:** the coastline is downloaded once from Natural Earth (via cartopy), so it
  needs internet access then.
- **Start-up:** about 10 s for the metadata and bin cells. The fine grid's profile
  positions (for the Profiles view) load in the background in ~30 s; a Profiles request
  before then waits for them once.
- **Speed after that:**

  | Change | Time |
  |---|---|
  | Recolor, change depth, change source | ~0.03–0.1 s |
  | Profiles view, whole record | ~1.3 s |
  | Profiles view, a few years | ~0.05 s |
  | Bin panels | ~0.1 s |
  | Sub-polygons, first view of an area | ~3 s (then cached) |

- **Memory:** about 3 GB.

### What the store holds

- **Per bin, per month and per source:** the mean T, S, Tclim and Sclim profiles on the
  97-level grid, and the number of profiles.
- **Per bin, all sources pooled and per source** (`<bin>/by_source/<source>/`):
  - the long-run mean and std of T − Tclim and S − Sclim at each depth;
  - for the sub-polygons, `T_anom_q` / `S_anom_q` (K × 97, float16) and `T_anom_n` /
    `S_anom_n` (valid count per depth). These hold every anomaly, sorted, if there are
    ≤ K at a depth, otherwise the (i + 0.5)/K quantiles. K is `--anom_quantiles`, default
    100.
- **Summary tables,** `summary/{fine,coarse}/{ALL,<source>}/`: every bin's anomaly mean /
  std and counts in one array each (`bin_ids`, `stats` (n × 4 × 97), `n` (n × 2 × 97)).
  The app loads a source's coloring from them in under a second.
- **SQLite profile table:** source, date, position, file and index of every profile, per
  bin.
- **Model equivalents** (`prof_Testim`, `prof_Sestim`), when a file has them (the chain
  doesn't produce them yet).

---

## Controls

| Control | Effect |
|---------|--------|
| Drag / scroll wheel | Pan / zoom (⌂ Reset view: back to the opening view) |
| Projection | World, Arctic or Southern Ocean |
| Show | Bins or Profiles |
| Grid resolution | Fine (10242) or coarse (2562) bins |
| Source | One source, or all combined |
| Color by (Bins) | Probability, or T / S anomaly mean / std at the chosen depth |
| Edges: anomaly std (Bins) | Std outlines on / off (anomaly-mean modes) |
| Probability (Bins) | A month counts with at least 1 / 2 / 3 profiles; out of all months or those since the bin's first profile (see Visit probability) |
| Time window (Profiles) | Years of the profiles drawn |
| Field | Temperature or salinity, for the bin panels |
| Anomaly reference | Subtract the climatology or the bin's long-run mean |
| Show source table | The bin's individual profiles, instead of the anomaly panel |

## Notes

- **Missing bins:** bins with no data for the selected source are not drawn.
- **Grid switch:** switching grid resolution clears the selected bin.
- **Coarse vs fine:** coarse bins show higher probabilities than fine bins, since each
  holds the profiles of about four fine bins.
- **Zoom limits:** the map can't zoom out past its opening view or in past 40×; the
  browser enforces this itself (`minscale` / `maxscale`). The polar maps pan freely
  (you can drag them well away from the pole); **⌂ Reset view** (top right of the map)
  returns to the opening zoom and orientation, keeping the depth, source and other
  settings.
- **Depth slider:** a tick for each of the 97 levels, depths written under every 16th.
  The chosen level and depth are shown above it; ◀ / ▶ step one level shallower / deeper.
- **Smooth zooming:** while you pan or zoom, the app sends the map nothing unless what it
  shows changes (new sub-polygons, the next outline-width step, a new set of dots).
- **Code layout:** `map_common.py` holds the projections, the bin-cell geometry, the
  coastline and the zoom / pan bookkeeping. `subpolygons.py` holds the sub-polygon and
  inset-outline geometry.
- **Older app variants** are in `x_old/` (not tracked by git).
