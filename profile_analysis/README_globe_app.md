# Profile globe apps

Interactive Dash/Plotly apps for exploring the NCEI-processed in-situ profiles
by instrument source, geodesic bin and time. They read a pre-built data store
(zarr + JSON + SQLite) made from one chain run by a build script.

| App | Shows | Store it reads |
|---|---|---|
| `app_allsource_globe.py` | Rotatable globe (or flat map) of geodesic bins coloured by sampling probability; click a bin for depth–time heatmaps of monthly-mean T / S (observed, climatology, anomaly) and a table of its profiles | `allsource_bin_timeseries_<ts>/`, or the anomaly store (the newest of either) |
| `app_allsource_anomaly_globe.py` | The same, plus colouring by the long-run anomaly mean / std (T or S, with a depth slider) | `allsource_bin_timeseries_with_anomalies_<ts>/` |
| `app_profile_scatter.py` | Individual profiles as dots coloured by source, with a date-range slider; click a dot for its bin's T / S timeseries | `allsource_bin_timeseries_with_anomalies_<ts>/` |

Sources are whatever the run holds, one per source folder of its `yearly/`
directory (today CTD_WOD, GLD_WOD, ITP, MEOP, MRB_WOD, PFL, XBT_WOD); files in
sub-folders (ITP/L2, ITP/L3) belong to their top folder's source. The
"All (combined)" view pools the sources, weighting each by its number of profiles.
Two grids: **fine** (10242 bins, `prof_bin_id_a`) and **coarse** (2562 bins,
`prof_bin_id_b`).

---

## Running an app

```bash
python profile_analysis/app_allsource_globe.py              # newest store
python profile_analysis/app_allsource_globe.py --data <store_dir>
```

Then open `http://127.0.0.1:8050`. All three apps use port 8050, so run one at a
time. Each picks the newest timestamped store in
`/Users/brucel/ecco/yip/profile_data/z_profile_file_analysis/` unless given
`--data`. Start-up loads only the JSON metadata; per-bin arrays are read from zarr,
and profile rows from SQLite, when you click a bin.

---

## Building a store from a chain run

```bash
R=/Users/brucel/ecco/yip/processed_by_NCEI_profile_data/profile_files_NCEI_processed_<ts>
python profile_analysis/build_allsource_bin_timeseries_with_anomalies.py --source_root $R
```

The anomaly store is a superset of the plain one, so this one build serves all
three apps (`build_allsource_bin_timeseries.py` still makes a plain store if you
want one; it is ~5 GB either way).

- `--source_root` can be the run folder: the scripts read `yearly/` inside it
  (older runs with source folders at the top level also work).
- Output: `z_profile_file_analysis/allsource_bin_timeseries[_with_anomalies]_<ts>/`,
  where `<ts>` is the run's timestamp, holding `<name>.zarr`,
  `<name>_metadata.json` and `<name>_profiles.db`. `--output <path>` overrides.
- The anomaly build can add its anomaly arrays to an existing store instead of
  rebuilding: `--append_anomalies`.
- Rebuild after every chain run you want to look at; an old store keeps showing
  the old run.
- If a build hangs after "Aggregating N bins...", rerun it with
  `OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1`.
- The geodesic bin-centre files are set at the top of the build scripts
  (`_DEFAULT_GEODESIC_FILE`, `_DEFAULT_GEODESIC_FILE_B`).

Per bin, per month and per source the store holds the mean T, S, Tclim, Sclim
profiles on the 97-level grid, the number of profiles, and (anomaly build) the
long-run mean and std of T − Tclim and S − Sclim at each depth. The model
equivalents (`prof_Testim`, `prof_Sestim`) are stored when a file has them; the
chain does not produce them yet.

---

## Controls

| Control | Effect |
|---------|--------|
| Drag globe | Rotate |
| Scroll wheel / zoom slider | Zoom |
| Grid resolution | Fine (10242) or coarse (2562) bins |
| Source | One source, or all combined |
| Probability threshold | Colour bins by P(≥1), P(≥2) or P(≥3) profiles per month |
| Probability basis | Denominator: the whole record, or since the bin's first profile |
| Color by (anomaly app) | Probability, or T / S anomaly mean / std at the chosen depth |
| Field | Temperature or salinity |
| Anomaly reference | Subtract the climatology or the bin's long-run mean |
| Click a bin | Depth–time heatmaps and (optionally) the bin's profile table |

## Notes

- Bins with no data for the selected source are not drawn.
- Switching grid resolution clears the selected bin.
- Coarse bins show higher probabilities than fine bins: each holds the profiles
  of about four fine bins.
- Known issues: rotation sticks near the poles (gimbal lock; a quaternion-based
  rotation would fix it); zooming out on the flat map is jerky.
- Older app variants are in `x_old/` (not tracked by git).
