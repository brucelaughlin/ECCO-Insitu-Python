"""
Preprocessing script for the multi-source interactive globe app.

Generalizes build_pfl_bin_timeseries.py to ingest EVERY post-NCEI source under
profile_files_NCEI_processed/ (CTD_WOD, GLD_WOD, MEOP, MRB_WOD, PFL, PFL_BGC,
XBT_WOD, ...), not just PFL. ITP is currently empty and is skipped with a log.

All sources share the same physical depth grid — each is a clean prefix of the
97-level reference grid (2, 4, 7, 10 ... m), differing only in max depth — so a
single pad-to-97 handles every source without interpolation.

Data is partitioned BY SOURCE so the app can either show one instrument in
isolation or pool them. The pooled ("all") view is derived on the fly; storing
per-source keeps every option open (the reverse is not recoverable).

Model equivalents (prof_Testim / prof_Sestim, the ECCO profiles-package
model-sampled T/S) are ingested WHEN PRESENT. As of this writing only MEOP has
them; the other sources have not had the upstream model-sampling step run yet.
The code fills estim for any source that has the fields and stores None-flagged
NaNs otherwise, so it auto-populates once upstream is run.

Per bin, per (year, month), per source, computes nanmean of:
  T, S           — observed profiles
  Tclim, Sclim   — climatology interpolated to each profile's location/time
  Testim, Sestim — model equivalent, where available

Pooling convention: the combined ("all") monthly mean is profile-weighted
(nanmean over every profile in the bin-month regardless of source), NOT a mean
of per-source means. A source contributing more profiles therefore weighs more.

STORAGE NOTES (why this pickle is compact):
  * Per-source arrays are the source of truth. The combined ("all") monthly
    arrays are NOT stored — the app shows one bin at a time and derives the
    profile-weighted combined view on click. Storing them for ~840k bin-months
    was pure redundancy.
  * Depth arrays are stored as float32 (ample for T/S) to halve array bytes.
  * Model equivalents (prof_Testim/Sestim) are currently ALL-NaN in every
    source file (the upstream ECCO model-sampling step has not populated them),
    so their arrays are NOT stored — only a per-source `has_estim` flag. Add the
    array storage back (see INGEST_ESTIM) once upstream fills the fields.

Output pickle (new file, does not overwrite the PFL one):
  'depth'            : np.ndarray float32 (97,)
  'sources'          : list[str]                sorted source names ingested
  'bins'             : dict[bin_id -> {
      'bin_id'       : int
      'months'       : dict[(year,month) -> {
          'n'        : int                      total profiles this month
          'by_source': dict[src -> {
              'T','S','Tclim','Sclim' : np.ndarray float32 (97,)
              'n'      : int
              'sources': list of {file, prof_idx, date, lon, lat, qual}
          }]
      }]
      'T_mean','S_mean' : np.ndarray float32 (97,)  combined long-run monthly mean
      'hit_counts'      : {1,2,3 -> int}         combined months with >= N profiles
      'first_ym','last_ym' : (year,month)
      'prob_1plus'      : float                  combined day-0 probability
      'by_source'       : dict[src -> {
          'T_mean','S_mean' : np.ndarray float32 (97,)
          'hit_counts'      : {1,2,3 -> int}
          'first_ym','last_ym' : (year,month)
          'n_profiles'      : int
          'has_estim'       : bool               True once upstream fills estim
      }]
  }]
  'bin_centres'      : dict[bin_id -> (lon, lat)]
  'dataset_first_ym' : (year, month)
  'dataset_last_ym'  : (year, month)
  'total_months'     : int
"""

import pickle
import numpy as np
import pandas as pd
import xarray as xr
from pathlib import Path
from collections import defaultdict

# ==============================================================================
# Paths
# ==============================================================================

SOURCE_ROOT   = Path('/Users/brucel/ecco/yip/profile_files_NCEI_processed')
GEODESIC_FILE = '/Users/brucel/ecco/yip/sample_data/ecco-insitu/sweet_gdrive/geodesic/10242_bin_locations.csv'
OUTPUT_PICKLE = Path('/Users/brucel/ecco/yip/profile_data/z_profile_file_analysis/allsource_bin_timeseries.pkl')

BIN_ID_VAR  = 'prof_bin_id_a'   # 10242-bin resolution, precomputed by NCEI chain
N_DEPTH_REF = 97                 # full reference depth grid length

# prof_Testim/Sestim are all-NaN in every source today (upstream model-sampling
# step not yet run). We still track a per-source has_estim flag. Flip this to
# True to also STORE the estim arrays once upstream populates them.
INGEST_ESTIM = False

# ==============================================================================
# Helpers
# ==============================================================================

def load_geodesic_centres(geodesic_file):
    df = pd.read_csv(geodesic_file, header=None, names=['lon', 'lat'])
    return {i + 1: (float(df['lon'].iloc[i]), float(df['lat'].iloc[i]))
            for i in range(len(df))}


def pad_to_ref(arr, n_ref=N_DEPTH_REF):
    """Pad a 1-D depth profile to n_ref levels with NaNs; stored as float32."""
    arr = np.asarray(arr, dtype=np.float32)
    if len(arr) == n_ref:
        return arr
    out = np.full(n_ref, np.nan, dtype=np.float32)
    out[:min(len(arr), n_ref)] = arr[:n_ref]
    return out


def stack_nanmean(profiles, key):
    """nanmean over a list of profile dicts for the given per-profile key."""
    return np.nanmean(np.stack([p[key] for p in profiles], axis=0),
                      axis=0).astype(np.float32)


def file_qualifier(stem):
    """Sub-source qualifier from the filename (e.g. PFL D/R/A, WOD OSD/noflag).

    The true source is the directory name; this is only a finer tag kept for
    the source table. Returns '' when there is no meaningful qualifier.
    """
    parts = stem.split('_')
    # parts like: ARGO_WO_1997_PFL_D__ncei_step_10  or  WOD_WO_1992_CTD_OSD__ncei
    # index 4 is the token after the instrument code, before the trailing __ncei
    if len(parts) > 4 and parts[4] and parts[4] != 'ncei':
        return parts[4]
    return ''

# ==============================================================================
# Reference depth grid
# ==============================================================================

print("Loading geodesic bin centres...")
bin_centres = load_geodesic_centres(GEODESIC_FILE)
print(f"  {len(bin_centres):,} bin centres loaded")

# Establish the 97-level reference depth grid from the first file that has it,
# searching across all sources.
DEPTH_REF = None
for _d in sorted(SOURCE_ROOT.iterdir()):
    if not _d.is_dir():
        continue
    for _f in sorted(_d.glob('*.nc')):
        try:
            _ds = xr.open_dataset(_f)
        except Exception:
            continue
        if _ds.sizes.get('iDEPTH', 0) == N_DEPTH_REF:
            DEPTH_REF = _ds['depth'].values.astype(float).ravel()
            _ds.close()
            break
        _ds.close()
    if DEPTH_REF is not None:
        break
if DEPTH_REF is None:
    raise RuntimeError(f"No post-NCEI file found with iDEPTH={N_DEPTH_REF}")
print(f"  Reference depth grid: {N_DEPTH_REF} levels, "
      f"{DEPTH_REF[0]:.0f}-{DEPTH_REF[-1]:.0f} m")

# ==============================================================================
# Ingest — bin_profiles[bin_id][source] = list of per-profile dicts
# ==============================================================================

source_dirs = sorted(d for d in SOURCE_ROOT.iterdir() if d.is_dir())
bin_profiles = defaultdict(lambda: defaultdict(list))
sources_seen = set()

for sdir in source_dirs:
    source = sdir.name
    files  = sorted(sdir.glob('*.nc'))
    if not files:
        print(f"\n=== {source} === no .nc files, skipping")
        continue
    print(f"\n=== {source} === {len(files)} files")

    for fi, fpath in enumerate(files):
        qual = file_qualifier(fpath.stem)
        try:
            ds = xr.open_dataset(fpath)
        except Exception as e:
            print(f"  [{fi+1}/{len(files)}] {fpath.name}: skipping ({e})")
            continue

        n_prof = ds.sizes.get('iPROF', 0)
        if n_prof == 0:
            ds.close()
            continue

        yyyymmdd = ds['prof_YYYYMMDD'].values.astype(int)
        years    = yyyymmdd // 10000
        months   = (yyyymmdd % 10000) // 100
        lons     = ds['prof_lon'].values.astype(float)
        lats     = ds['prof_lat'].values.astype(float)
        T_all    = ds['prof_T'].values.astype(float)
        S_all    = ds['prof_S'].values.astype(float)
        Tc_all   = ds['prof_Tclim'].values.astype(float)
        Sc_all   = ds['prof_Sclim'].values.astype(float)
        bin_ids  = ds[BIN_ID_VAR].values.astype(float)

        estim_present = 'prof_Testim' in ds.variables and 'prof_Sestim' in ds.variables
        if estim_present:
            Te_all = ds['prof_Testim'].values.astype(float)
            Se_all = ds['prof_Sestim'].values.astype(float)

        for i in range(n_prof):
            bid = bin_ids[i]
            if np.isnan(bid) or bid <= 0:
                continue
            # estim is "real" only if the field exists AND is not all-NaN for
            # this profile (today it is all-NaN everywhere; see INGEST_ESTIM)
            rec_has_estim = estim_present and not np.all(np.isnan(Te_all[i]))
            rec = {
                'year':     int(years[i]),
                'month':    int(months[i]),
                'file':     fpath.name,
                'prof_idx': i,
                'date':     int(yyyymmdd[i]),
                'lon':      float(lons[i]),
                'lat':      float(lats[i]),
                'qual':     qual,
                'T':        pad_to_ref(T_all[i]),
                'S':        pad_to_ref(S_all[i]),
                'Tclim':    pad_to_ref(Tc_all[i]),
                'Sclim':    pad_to_ref(Sc_all[i]),
                'has_estim': rec_has_estim,
            }
            if INGEST_ESTIM and rec_has_estim:
                rec['Testim'] = pad_to_ref(Te_all[i])
                rec['Sestim'] = pad_to_ref(Se_all[i])
            bin_profiles[int(bid)][source].append(rec)

        sources_seen.add(source)
        ds.close()
        if (fi + 1) % 20 == 0 or fi + 1 == len(files):
            print(f"  [{fi+1}/{len(files)}] ingested; {len(bin_profiles):,} bins so far",
                  flush=True)

sources_sorted = sorted(sources_seen)
print(f"\nSources ingested: {', '.join(sources_sorted)}")
print(f"Aggregating {len(bin_profiles):,} bins...")

# ==============================================================================
# Dataset-wide calendar span (probability denominator)
# ==============================================================================

all_ym = set()
for src_map in bin_profiles.values():
    for profs in src_map.values():
        for p in profs:
            all_ym.add((p['year'], p['month']))
all_ym = sorted(all_ym)
if all_ym:
    (y0, m0), (y1, m1) = all_ym[0], all_ym[-1]
    total_months = (y1 - y0) * 12 + (m1 - m0) + 1
else:
    y0 = m0 = y1 = m1 = None
    total_months = 1


def hit_counts_from_groups(ym_groups):
    """{1,2,3 -> #months with >= N profiles} from a (y,m)->list mapping."""
    hc = {1: 0, 2: 0, 3: 0}
    for group in ym_groups.values():
        n = len(group)
        if n >= 1: hc[1] += 1
        if n >= 2: hc[2] += 1
        if n >= 3: hc[3] += 1
    return hc

# ==============================================================================
# Aggregate per bin
# ==============================================================================

bins_out = {}

for bid, src_map in bin_profiles.items():
    # ---- combined (all sources pooled, profile-weighted) ----
    combined_ym = defaultdict(list)
    for profs in src_map.values():
        for p in profs:
            combined_ym[(p['year'], p['month'])].append(p)

    months_dict   = {}
    all_T_monthly = []
    all_S_monthly = []

    for ym in sorted(combined_ym):
        group = combined_ym[ym]
        # combined monthly mean tracked only to build the long-run T_mean/S_mean;
        # NOT stored per-month (app derives the combined view on click)
        all_T_monthly.append(stack_nanmean(group, 'T'))
        all_S_monthly.append(stack_nanmean(group, 'S'))

        # per-source breakdown for this month, built directly from src_map so
        # source identity stays exact (combined `group` has lost its tags)
        by_source = {}
        for source, profs in src_map.items():
            grp = [p for p in profs if (p['year'], p['month']) == ym]
            if not grp:
                continue
            entry = {
                'T':      stack_nanmean(grp, 'T'),
                'S':      stack_nanmean(grp, 'S'),
                'Tclim':  stack_nanmean(grp, 'Tclim'),
                'Sclim':  stack_nanmean(grp, 'Sclim'),
                'n':      len(grp),
                'sources': [{
                    'file':     p['file'],
                    'prof_idx': p['prof_idx'],
                    'date':     p['date'],
                    'lon':      p['lon'],
                    'lat':      p['lat'],
                    'qual':     p['qual'],
                } for p in grp],
            }
            if INGEST_ESTIM and any(p['has_estim'] for p in grp):
                entry['Testim'] = stack_nanmean([p for p in grp if p['has_estim']], 'Testim')
                entry['Sestim'] = stack_nanmean([p for p in grp if p['has_estim']], 'Sestim')
            by_source[source] = entry

        months_dict[ym] = {
            'n':         len(group),
            'by_source': by_source,
        }

    T_long_mean = (np.nanmean(np.stack(all_T_monthly, axis=0), axis=0).astype(np.float32)
                   if all_T_monthly else np.full(N_DEPTH_REF, np.nan, dtype=np.float32))
    S_long_mean = (np.nanmean(np.stack(all_S_monthly, axis=0), axis=0).astype(np.float32)
                   if all_S_monthly else np.full(N_DEPTH_REF, np.nan, dtype=np.float32))

    combined_hits = hit_counts_from_groups(combined_ym)
    bin_yms       = sorted(combined_ym.keys())

    # ---- per-source rollups ----
    by_source_roll = {}
    for source, profs in src_map.items():
        src_ym = defaultdict(list)
        for p in profs:
            src_ym[(p['year'], p['month'])].append(p)
        src_T_monthly = [stack_nanmean(g, 'T') for g in
                         (src_ym[k] for k in sorted(src_ym))]
        src_S_monthly = [stack_nanmean(g, 'S') for g in
                         (src_ym[k] for k in sorted(src_ym))]
        src_yms = sorted(src_ym.keys())
        by_source_roll[source] = {
            'T_mean':     (np.nanmean(np.stack(src_T_monthly, axis=0), axis=0).astype(np.float32)
                           if src_T_monthly else np.full(N_DEPTH_REF, np.nan, dtype=np.float32)),
            'S_mean':     (np.nanmean(np.stack(src_S_monthly, axis=0), axis=0).astype(np.float32)
                           if src_S_monthly else np.full(N_DEPTH_REF, np.nan, dtype=np.float32)),
            'hit_counts': hit_counts_from_groups(src_ym),
            'first_ym':   src_yms[0],
            'last_ym':    src_yms[-1],
            'n_profiles': len(profs),
            'has_estim':  any(p['has_estim'] for p in profs),
        }

    bins_out[bid] = {
        'bin_id':     bid,
        'months':     months_dict,
        'T_mean':     T_long_mean,
        'S_mean':     S_long_mean,
        'hit_counts': combined_hits,
        'first_ym':   bin_yms[0],
        'last_ym':    bin_yms[-1],
        'prob_1plus': combined_hits[1] / total_months,
        'by_source':  by_source_roll,
    }

# ==============================================================================
# Save
# ==============================================================================

result = {
    'depth':            DEPTH_REF.astype(np.float32),
    'sources':          sources_sorted,
    'bins':             bins_out,
    'bin_centres':      bin_centres,
    'dataset_first_ym': (y0, m0) if all_ym else None,
    'dataset_last_ym':  (y1, m1) if all_ym else None,
    'total_months':     total_months,
}

print(f"\nSaving to {OUTPUT_PICKLE}...")
OUTPUT_PICKLE.parent.mkdir(parents=True, exist_ok=True)
with open(OUTPUT_PICKLE, 'wb') as f:
    pickle.dump(result, f, protocol=4)

size_mb = OUTPUT_PICKLE.stat().st_size / 1e6
print(f"Done. {len(bins_out):,} bins across {len(sources_sorted)} sources. "
      f"File size: {size_mb:.1f} MB")
