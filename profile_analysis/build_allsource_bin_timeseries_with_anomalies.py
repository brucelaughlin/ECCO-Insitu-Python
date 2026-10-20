"""
Builds the data store for the profile map app (app_profile_map.py) from one NCEI run.

Per bin, month and source: mean T / S / Tclim / Sclim profiles; per bin, per-depth
anomaly statistics (mean and std of T - Tclim and S - Sclim across all profiles ever
in the bin) and the anomaly quantiles for the sub-polygons.
These are stored directly in the zarr store alongside T_mean/S_mean:

  fine/<bid>/T_anom_mean   (97,) float32
  fine/<bid>/T_anom_std    (97,) float32
  fine/<bid>/S_anom_mean   (97,) float32
  fine/<bid>/S_anom_std    (97,) float32
  fine/<bid>/T_anom_q      (K, 97) float16   anomalies at each depth: all of them, sorted, if <= K,
  fine/<bid>/S_anom_q      (K, 97) float16   else the (i+0.5)/K quantiles (K = --anom_quantiles, 100);
  fine/<bid>/T_anom_n      (97,) int32       for the map app's sub-polygons
  fine/<bid>/S_anom_n      (97,) int32       valid anomaly count per depth
  fine/<bid>/by_source/<src>/{T,S}_anom_{mean,std,q,n}   the same, from that source's profiles only
  summary/{fine,coarse}/{ALL,<src>}/{bin_ids (n,), stats (n, 4, 97), n (n, 2, 97)}
                                             every bin's anomaly mean/std (T, S) and counts in one
                                             table per source, so the map app loads them at once

Run modes
---------
Default (full rebuild):
    python build_allsource_bin_timeseries_with_anomalies.py --source_root <dir>

Append anomaly arrays to an existing store (skips ingest + aggregate):
    python build_allsource_bin_timeseries_with_anomalies.py \\
        --source_root <dir> --append_anomalies

Output
------
  <output>.zarr              zarr store (arrays)
  <output>_metadata.json     JSON (scalars, profile records, bin centers)

Default output path:
  /Users/brucel/ecco/yip/profile_data/z_profile_file_analysis/
      allsource_bin_timeseries_with_anomalies
"""

import argparse
import sys
import concurrent.futures
import json
import multiprocessing as mp
import os
os.environ.setdefault('OMP_NUM_THREADS', '1')
os.environ.setdefault('OPENBLAS_NUM_THREADS', '1')
import re
import sqlite3
import warnings
import numpy as np
import pandas as pd
import xarray as xr
import zarr
from pathlib import Path
from collections import defaultdict

# ==============================================================================
# Paths
# ==============================================================================

# Defaults for this project's machine; other users pass --geodesic_dir and --output.
_DEFAULT_GEODESIC_DIR = '/Users/brucel/ecco/yip/sample_data/ecco-insitu/sweet_gdrive/geodesic'
_GEODESIC_NAMES       = ('10242_bin_locations.csv', '02562_bin_locations.csv')   # fine, coarse
_OUTPUT_BASE          = ('/Users/brucel/ecco/yip/profile_data/z_profile_file_analysis/'
                         'allsource_bin_timeseries_with_anomalies')

# nanmean/nanstd over all-NaN depth levels (e.g. XBT shallower than 97 levels) is expected.
warnings.filterwarnings('ignore', message='Mean of empty slice', category=RuntimeWarning)
warnings.filterwarnings('ignore', message='Degrees of freedom <= 0', category=RuntimeWarning)

_parser = argparse.ArgumentParser(
    description='Build allsource bin timeseries with anomaly stats (zarr + JSON).')
_parser.add_argument('--source_root', required=True,
                     help='Post-NCEI profile files: a chain run folder (its yearly/ is used), or any '
                          'folder with one sub-folder per source holding the .nc files')
_parser.add_argument('--geodesic_dir', default=_DEFAULT_GEODESIC_DIR,
                     help='Folder with the geodesic bin centers, ' + ' and '.join(_GEODESIC_NAMES)
                          + ' (ECCO in-situ support data)')
_parser.add_argument('--output', default=None,
                     help='Store folder to write (it holds <name>.zarr, <name>_metadata.json and '
                          '<name>_profiles.db). Default: allsource_bin_timeseries_with_anomalies_<ts>, '
                          '<ts> from the run folder name, in this project\'s store folder if it exists, '
                          'else in the current folder.')
_parser.add_argument('--append_anomalies', action='store_true',
                     help='Skip ingest/aggregate; open existing store and write only the '
                          'T/S_anom_mean/std (and quantile) arrays. Requires the store to already exist.')
_parser.add_argument('--anom_quantiles', type=int, default=100,
                     help='K: per bin and depth, store the anomalies (T-Tclim, S-Sclim) as at most K '
                          'values for the map app\'s sub-polygons: all of them, sorted, when there are '
                          '<= K, else the (i+0.5)/K quantiles (default 100)')
_args = _parser.parse_args()
ANOM_Q = _args.anom_quantiles

SOURCE_ROOT      = Path(_args.source_root)
GEODESIC_FILE, GEODESIC_FILE_B = (str(Path(_args.geodesic_dir) / n) for n in _GEODESIC_NAMES)
APPEND_ANOMALIES = _args.append_anomalies
for _f in (GEODESIC_FILE, GEODESIC_FILE_B):
    if not Path(_f).is_file():
        sys.exit(f'Geodesic bin-center file not found: {_f}\n(pass --geodesic_dir <folder with '
                 + ' and '.join(_GEODESIC_NAMES) + '>)')
if not SOURCE_ROOT.is_dir():
    sys.exit(f'--source_root not found: {SOURCE_ROOT}')
if not Path(_OUTPUT_BASE).parent.is_dir():          # not this project's machine
    _OUTPUT_BASE = str(Path.cwd() / Path(_OUTPUT_BASE).name)

def _derive_output(source_root, base):
    """Append _YYYYMMDD_HHMMSS timestamp from source_root name to base, if present."""
    m = re.search(r'_(\d{8}_\d{6})$', source_root.name)
    if m:
        return base + '_' + m.group(1)
    return base

OUTPUT_DIR = Path(_args.output) if _args.output else Path(_derive_output(SOURCE_ROOT, _OUTPUT_BASE))
print(f"Output path stem: {OUTPUT_DIR}")

# Run folders hold yearly/ (per-year files, one sub-folder per source) and
# ecco_arrays/; read the sources from yearly/.  (Older runs: sources at the top.)
if (SOURCE_ROOT / 'yearly').is_dir():
    SOURCE_ROOT = SOURCE_ROOT / 'yearly'


def _source_dirs(root):
    """[(source_name, folder)]: one source per top-level folder holding .nc files,
    directly or in sub-folders (ITP/L2 and ITP/L3 together form source ITP)."""
    return [(d.name, d) for d in sorted(p for p in root.iterdir() if p.is_dir()) if any(d.rglob('*.nc'))]
    print(f"Reading sources from {SOURCE_ROOT}")

BIN_ID_VAR   = 'prof_bin_id_a'   # 10242-bin resolution
BIN_ID_VAR_B = 'prof_bin_id_b'   # 2562-bin resolution
N_DEPTH_REF  = 97

INGEST_ESTIM = False

_INGEST_WORKERS = min(8, max(1, (os.cpu_count() or 1) - 1))

# Output layout: a directory named after the stem, containing both files.
#   <stem>/
#       <stem>.zarr
#       <stem>_metadata.json
OUTPUT_CONTAINER = OUTPUT_DIR
OUTPUT_ZARR      = OUTPUT_CONTAINER / (OUTPUT_DIR.name + '.zarr')
OUTPUT_JSON      = OUTPUT_CONTAINER / (OUTPUT_DIR.name + '_metadata.json')

# ==============================================================================
# Helpers
# ==============================================================================

def load_geodesic_centers(geodesic_file):
    df = pd.read_csv(geodesic_file, header=None, names=['lon', 'lat'])
    return {i + 1: (float(df['lon'].iloc[i]), float(df['lat'].iloc[i]))
            for i in range(len(df))}


def pad_to_ref(arr, n_ref=N_DEPTH_REF):
    arr = np.asarray(arr, dtype=np.float32)
    if len(arr) == n_ref:
        return arr
    out = np.full(n_ref, np.nan, dtype=np.float32)
    out[:min(len(arr), n_ref)] = arr[:n_ref]
    return out


def stack_nanmean(profiles, key):
    return np.nanmean(np.stack([p[key] for p in profiles], axis=0),
                      axis=0).astype(np.float32)


def file_qualifier(stem):
    parts = stem.split('_')
    if len(parts) > 4 and parts[4] and parts[4] != 'ncei':
        return parts[4]
    return ''


def _ym_str(ym):
    return f"{ym[0]}-{ym[1]:02d}"


def _ym_int(ym):
    return ym[0] * 100 + ym[1]

# ==============================================================================
# Anomaly-only append mode
# ==============================================================================

def _anom_quantiles(stack, k=None):
    """Compact form of a bin's anomalies (n_prof, N_DEPTH) for the map app's sub-polygons.

    Returns q (k, N_DEPTH) float16 and n (N_DEPTH,) int32, the valid count per depth.
    At each depth, q holds the sorted valid values when there are <= k of them (the
    rest NaN), else the (i + 0.5) / k quantiles, i = 0..k-1.
    """
    k = ANOM_Q if k is None else k
    n = np.isfinite(stack).sum(axis=0).astype(np.int32)
    q = np.full((k, stack.shape[1]), np.nan, dtype=np.float32)
    srt = np.sort(stack, axis=0)                                # NaNs sort last
    probs = (np.arange(k) + 0.5) / k
    for d in np.flatnonzero(n):
        v = srt[:n[d], d]
        if n[d] <= k:
            q[:n[d], d] = v
        else:
            q[:, d] = np.quantile(v, probs)
    return q.astype(np.float16), n


def _anom_q_arrays(T_anom, S_anom):
    Tq, Tn = _anom_quantiles(T_anom)
    Sq, Sn = _anom_quantiles(S_anom)
    return {'T_anom_q': Tq, 'T_anom_n': Tn, 'S_anom_q': Sq, 'S_anom_n': Sn}


_ANOM_Q_KEYS = ('T_anom_q', 'T_anom_n', 'S_anom_q', 'S_anom_n')


ANOM_ALL = 'ALL'                 # key for "all sources pooled" in the per-source anomaly dicts


def _anom_arrays(profs):
    """Anomaly arrays (T-Tclim, S-Sclim) of a set of profiles: mean, std, quantiles, counts."""
    T = np.stack([p['T'] - p['Tclim'] for p in profs], axis=0).astype(np.float32)
    S = np.stack([p['S'] - p['Sclim'] for p in profs], axis=0).astype(np.float32)
    return {'T_anom_mean': np.nanmean(T, axis=0).astype(np.float32),
            'T_anom_std':  np.nanstd( T, axis=0).astype(np.float32),
            'S_anom_mean': np.nanmean(S, axis=0).astype(np.float32),
            'S_anom_std':  np.nanstd( S, axis=0).astype(np.float32),
            **_anom_q_arrays(T, S)}


def _anom_by_source(src_map):
    """{ANOM_ALL: arrays of all the bin's profiles, src: arrays of that source's profiles}."""
    out = {ANOM_ALL: _anom_arrays([p for profs in src_map.values() for p in profs])}
    for src, profs in src_map.items():
        if profs:
            out[src] = _anom_arrays(profs)
    return out


def _write_anom(group, arrays):
    for name, arr in arrays.items():
        if name in group:
            del group[name]
        group.create_array(name, data=arr, chunks=arr.shape)


def _summary_rows(anom):
    """The small per-bin pieces kept for the summary tables: stats (4, D), counts (2, D)."""
    return (np.stack([anom['T_anom_mean'], anom['T_anom_std'], anom['S_anom_mean'], anom['S_anom_std']]),
            np.stack([anom['T_anom_n'], anom['S_anom_n']]).astype(np.int32))


def _write_summary(store, resolution_label, rows_by_bin):
    """summary/<res>/<ALL | source>/{bin_ids (n,), stats (n, 4, D), n (n, 2, D)}: one read
    gives the map app every bin's anomaly mean / std / counts for that source.
    rows_by_bin: {bid: {key: (stats, counts)}}"""
    sg = store.require_group('summary').require_group(resolution_label)
    keys = sorted({k for rows in rows_by_bin.values() for k in rows})
    for key in keys:
        bids = sorted(b for b, rows in rows_by_bin.items() if key in rows)
        _write_anom(sg.require_group(key), {
            'bin_ids': np.array(bids, dtype=np.int32),
            'stats':   np.stack([rows_by_bin[b][key][0] for b in bids]).astype(np.float32),
            'n':       np.stack([rows_by_bin[b][key][1] for b in bids]).astype(np.int32)})
    print(f"  summary/{resolution_label}: {len(keys)} tables ({', '.join(keys)})", flush=True)


def _append_anomaly_arrays_from_store(store, resolution_label, bin_profiles_res):
    """Compute and write the anomaly arrays (mean, std, quantiles, counts; pooled at the
    bin root and per source under by_source/<src>) and the summary tables into an
    already-open zarr store."""
    grp  = store[resolution_label]
    bids = [b for b in bin_profiles_res if any(bin_profiles_res[b].values())]
    print(f"  Writing anomaly arrays for {len(bids):,} {resolution_label} bins...", flush=True)

    def one(bid):
        anom = _anom_by_source(bin_profiles_res[bid])
        bg = grp.require_group(str(bid))
        _write_anom(bg, anom[ANOM_ALL])
        for src, arrays in anom.items():
            if src != ANOM_ALL:
                _write_anom(bg.require_group('by_source').require_group(src), arrays)
        return bid, {k: _summary_rows(a) for k, a in anom.items()}

    rows_by_bin = {}
    with concurrent.futures.ThreadPoolExecutor(max_workers=max(1, (os.cpu_count() or 1) - 1)) as pool:
        for i, (bid, rows) in enumerate(pool.map(one, bids)):
            rows_by_bin[bid] = rows
            if i % 1000 == 0:
                print(f"    {i}/{len(bids)}", flush=True)
    _write_summary(store, resolution_label, rows_by_bin)

# ==============================================================================
# Ingest worker
# ==============================================================================

def _ingest_file(fpath_source):
    fpath, source = fpath_source
    qual = file_qualifier(fpath.stem)
    try:
        ds = xr.open_dataset(fpath)
    except Exception as e:
        return source, fpath.name, None, None, str(e)

    n_prof = ds.sizes.get('iPROF', 0)
    if n_prof == 0:
        ds.close()
        return source, fpath.name, [], [], None

    yyyymmdd = ds['prof_YYYYMMDD'].values.astype(int)
    years    = yyyymmdd // 10000
    months   = (yyyymmdd % 10000) // 100
    lons     = ds['prof_lon'].values.astype(float)
    lats     = ds['prof_lat'].values.astype(float)
    n_depth  = ds.sizes.get('iDEPTH', N_DEPTH_REF)
    _nan2d   = np.full((n_prof, n_depth), np.nan, dtype=float)
    T_all    = ds['prof_T'].values.astype(float)      if 'prof_T'     in ds else _nan2d
    S_all    = ds['prof_S'].values.astype(float)      if 'prof_S'     in ds else _nan2d
    Tc_all   = ds['prof_Tclim'].values.astype(float)  if 'prof_Tclim' in ds else _nan2d
    Sc_all   = ds['prof_Sclim'].values.astype(float)  if 'prof_Sclim' in ds else _nan2d
    bin_ids   = ds[BIN_ID_VAR].values.astype(float)
    bin_ids_b = ds[BIN_ID_VAR_B].values.astype(float) if BIN_ID_VAR_B in ds else None

    estim_present = 'prof_Testim' in ds.variables and 'prof_Sestim' in ds.variables
    if estim_present:
        Te_all = ds['prof_Testim'].values.astype(float)
        Se_all = ds['prof_Sestim'].values.astype(float)

    recs_fine   = []
    recs_coarse = []
    for i in range(n_prof):
        bid = bin_ids[i]
        if np.isnan(bid) or bid <= 0:
            continue
        rec_has_estim = estim_present and not np.all(np.isnan(Te_all[i]))
        rec = {
            'year':     int(years[i]),
            'month':    int(months[i]),
            'file':     str(fpath.relative_to(SOURCE_ROOT / source)),   # L2/..., L3/... for ITP
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
        recs_fine.append((int(bid), rec))

        if bin_ids_b is not None:
            bid_b = bin_ids_b[i]
            if not (np.isnan(bid_b) or bid_b <= 0):
                recs_coarse.append((int(bid_b), rec))

    ds.close()
    return source, fpath.name, recs_fine, recs_coarse, None

# ==============================================================================
# Aggregate
# ==============================================================================

def hit_counts_from_groups(ym_groups):
    hc = {1: 0, 2: 0, 3: 0}
    for group in ym_groups.values():
        n = len(group)
        if n >= 1: hc[1] += 1
        if n >= 2: hc[2] += 1
        if n >= 3: hc[3] += 1
    return hc


def _agg_bin(bid_src_map_tm):
    """Aggregate a single bin. Top-level so ProcessPoolExecutor can pickle it."""
    bid, src_map, tm = bid_src_map_tm

    combined_ym = defaultdict(list)
    for profs in src_map.values():
        for p in profs:
            combined_ym[(p['year'], p['month'])].append(p)

    src_ym_map = {}
    for source, profs in src_map.items():
        grp_by_ym = defaultdict(list)
        for p in profs:
            grp_by_ym[(p['year'], p['month'])].append(p)
        src_ym_map[source] = grp_by_ym

    months_dict   = {}
    all_T_monthly = []
    all_S_monthly = []

    for ym in sorted(combined_ym):
        group = combined_ym[ym]
        all_T_monthly.append(stack_nanmean(group, 'T'))
        all_S_monthly.append(stack_nanmean(group, 'S'))

        by_source = {}
        for source in src_map:
            grp = src_ym_map[source].get(ym, [])
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

        months_dict[ym] = {'n': len(group), 'by_source': by_source}

    T_long_mean = (np.nanmean(np.stack(all_T_monthly, axis=0), axis=0).astype(np.float32)
                   if all_T_monthly else np.full(N_DEPTH_REF, np.nan, dtype=np.float32))
    S_long_mean = (np.nanmean(np.stack(all_S_monthly, axis=0), axis=0).astype(np.float32)
                   if all_S_monthly else np.full(N_DEPTH_REF, np.nan, dtype=np.float32))

    anom = _anom_by_source(src_map)          # pooled (ANOM_ALL) and per source

    combined_hits = hit_counts_from_groups(combined_ym)
    bin_yms       = sorted(combined_ym.keys())

    by_source_roll = {}
    for source, profs in src_map.items():
        src_ym  = src_ym_map[source]
        src_yms = sorted(src_ym.keys())
        src_T_monthly = [stack_nanmean(src_ym[k], 'T') for k in src_yms]
        src_S_monthly = [stack_nanmean(src_ym[k], 'S') for k in src_yms]
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

    return bid, {
        'bin_id':      bid,
        'months':      months_dict,
        'T_mean':      T_long_mean,
        'S_mean':      S_long_mean,
        **anom[ANOM_ALL],
        'anom_by_source': {k: v for k, v in anom.items() if k != ANOM_ALL},
        'hit_counts':  combined_hits,
        'first_ym':    bin_yms[0],
        'last_ym':     bin_yms[-1],
        'prob_1plus':  combined_hits[1] / tm,
        'by_source':   by_source_roll,
    }


_AGG_WORKERS = max(1, (os.cpu_count() or 1) - 1)


def aggregate_bins(bp, tm):
    """Aggregate bin_profiles in parallel across CPUs (bins are independent)."""
    work = [(bid, dict(src_map), tm) for bid, src_map in bp.items()]
    out  = {}
    done = 0
    # 'fork' copies the parent process (with bin_profiles already loaded) into
    # each worker instantly — avoiding the cost of re-importing and re-ingesting
    # that 'spawn' (the macOS default) would cause.
    #
    # KNOWN RISK: fork + OpenMP/BLAS can deadlock if a numpy thread pool was
    # active at fork time (the child inherits a locked mutex). Symptoms: workers
    # hang silently, no output after "Aggregating N bins...".
    # FIX: run with OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 python build_allsource...
    # or switch mp_context to 'spawn' and accept the re-import overhead.
    ctx  = mp.get_context('fork')
    with concurrent.futures.ProcessPoolExecutor(max_workers=_AGG_WORKERS,
                                                mp_context=ctx) as pool:
        futures = {pool.submit(_agg_bin, item): item for item in work}
        for fut in concurrent.futures.as_completed(futures):
            bid, bd = fut.result()
            out[bid] = bd
            done += 1
            if done % 1000 == 0:
                print(f"  aggregated {done}/{len(work)}", flush=True)
    return out

# ==============================================================================
# Zarr + JSON serialisation
# ==============================================================================

def write_bin_arrays(parent_group, bid, bd):
    bg = parent_group.require_group(str(bid))
    for key in ('T_mean', 'S_mean', 'T_anom_mean', 'T_anom_std', 'S_anom_mean', 'S_anom_std'):
        bg.create_array(key, data=np.asarray(bd[key], dtype=np.float32), chunks=(N_DEPTH_REF,))
    for key in _ANOM_Q_KEYS:                       # (K, N_DEPTH) float16 / (N_DEPTH,) int32
        bg.create_array(key, data=bd[key], chunks=bd[key].shape)

    sorted_yms = sorted(bd['months'].keys())
    nm         = len(sorted_yms)
    ym_ints    = np.array([_ym_int(ym) for ym in sorted_yms], dtype=np.int32)
    bg.create_array('ym_index', data=ym_ints, chunks=(nm,))

    nan_row = np.full(N_DEPTH_REF, np.nan, dtype=np.float32)
    bsg = bg.require_group('by_source')
    for src, roll in bd['by_source'].items():
        rg = bsg.require_group(src)
        rg.create_array('T_mean', data=np.asarray(roll['T_mean'], dtype=np.float32), chunks=(N_DEPTH_REF,))
        rg.create_array('S_mean', data=np.asarray(roll['S_mean'], dtype=np.float32), chunks=(N_DEPTH_REF,))
        if src in bd.get('anom_by_source', {}):
            _write_anom(rg, bd['anom_by_source'][src])      # this source's anomaly mean/std/quantiles
        for key in ('T', 'S', 'Tclim', 'Sclim'):
            mat = np.stack([
                np.asarray(bd['months'][ym]['by_source'][src][key], dtype=np.float32)
                if (src in bd['months'][ym]['by_source']) else nan_row
                for ym in sorted_yms
            ], axis=0)   # (n_months, N_DEPTH_REF)
            rg.create_array(key, data=mat, chunks=(nm, N_DEPTH_REF))


def bin_meta(bd):
    months_meta = {}
    for ym, m in bd['months'].items():
        by_src = {}
        for src, e in m['by_source'].items():
                by_src[src] = {'n': e['n']}   # profile records moved to SQLite
        months_meta[_ym_str(ym)] = {'n': m['n'], 'by_source': by_src}

    by_src_roll = {}
    for src, roll in bd['by_source'].items():
        by_src_roll[src] = {
            'hit_counts': {str(k): v for k, v in roll['hit_counts'].items()},
            'first_ym':   list(roll['first_ym']),
            'last_ym':    list(roll['last_ym']),
            'n_profiles': roll['n_profiles'],
            'has_estim':  roll['has_estim'],
        }

    return {
        'bin_id':     bd['bin_id'],
        'hit_counts': {str(k): v for k, v in bd['hit_counts'].items()},
        'first_ym':   list(bd['first_ym']),
        'last_ym':    list(bd['last_ym']),
        'prob_1plus': float(bd['prob_1plus']),
        'months':     months_meta,
        'by_source':  by_src_roll,
    }


def write_zarr_and_json(bins_fine, bins_coarse_dict, bin_centers_fine,
                        bin_centers_coarse_dict, depth_arr, sources,
                        dataset_first_ym, dataset_last_ym, total_months):

    print(f"\nWriting zarr store → {OUTPUT_ZARR} ...")
    store = zarr.open_group(str(OUTPUT_ZARR), mode='w')
    store.create_array('depth', data=depth_arr.astype(np.float32), chunks=(N_DEPTH_REF,))

    fine_g   = store.require_group('fine')
    coarse_g = store.require_group('coarse')

    _WRITE_WORKERS = max(1, (os.cpu_count() or 1) - 1)

    def _write_one(args):
        parent_group, bid, bd = args
        write_bin_arrays(parent_group, bid, bd)
        return bid

    n_fine = len(bins_fine)
    done   = 0
    with concurrent.futures.ThreadPoolExecutor(max_workers=_WRITE_WORKERS) as pool:
        futs = {pool.submit(_write_one, (fine_g, bid, bd)): bid
                for bid, bd in bins_fine.items()}
        for fut in concurrent.futures.as_completed(futs):
            fut.result()
            done += 1
            if done % 1000 == 0:
                print(f"  fine bins: {done}/{n_fine}", flush=True)
    print(f"  fine bins: {n_fine}/{n_fine}", flush=True)

    n_coarse = len(bins_coarse_dict)
    done     = 0
    with concurrent.futures.ThreadPoolExecutor(max_workers=_WRITE_WORKERS) as pool:
        futs = {pool.submit(_write_one, (coarse_g, bid, bd)): bid
                for bid, bd in bins_coarse_dict.items()}
        for fut in concurrent.futures.as_completed(futs):
            fut.result()
            done += 1
            if done % 500 == 0:
                print(f"  coarse bins: {done}/{n_coarse}", flush=True)
    print(f"  coarse bins: {n_coarse}/{n_coarse}", flush=True)

    for label, bdict in (('fine', bins_fine), ('coarse', bins_coarse_dict)):
        _write_summary(store, label, {
            bid: {ANOM_ALL: _summary_rows(bd), **{s: _summary_rows(a) for s, a in bd.get('anom_by_source', {}).items()}}
            for bid, bd in bdict.items()})

    zarr_size_mb = sum(f.stat().st_size for f in OUTPUT_ZARR.rglob('*') if f.is_file()) / 1e6
    print(f"  zarr done — {zarr_size_mb:.1f} MB")

    print(f"Writing JSON metadata → {OUTPUT_JSON} ...")
    meta = {
        'depth_list':         depth_arr.astype(np.float32).tolist(),
        'sources':            sources,
        'dataset_first_ym':   list(dataset_first_ym) if dataset_first_ym else None,
        'dataset_last_ym':    list(dataset_last_ym)  if dataset_last_ym  else None,
        'total_months':       total_months,
        'bin_centers':        {str(k): list(v) for k, v in bin_centers_fine.items()},
        'bin_centers_coarse': {str(k): list(v) for k, v in bin_centers_coarse_dict.items()},
        'bins_meta':          {str(bid): bin_meta(bd) for bid, bd in bins_fine.items()},
        'bins_coarse_meta':   {str(bid): bin_meta(bd) for bid, bd in bins_coarse_dict.items()},
    }
    with open(OUTPUT_JSON, 'w') as f:
        json.dump(meta, f, separators=(',', ':'))
    json_size_mb = OUTPUT_JSON.stat().st_size / 1e6
    print(f"  JSON done — {json_size_mb:.1f} MB")

# ==============================================================================
# Main
# ==============================================================================

print("Loading geodesic bin centers...")
bin_centers   = load_geodesic_centers(GEODESIC_FILE)
bin_centers_b = load_geodesic_centers(GEODESIC_FILE_B)
print(f"  {len(bin_centers):,} fine (10242) + {len(bin_centers_b):,} coarse (2562) centers loaded")

# --append_anomalies: open existing store, reuse existing JSON, write only anomaly arrays.
if APPEND_ANOMALIES:
    if not OUTPUT_ZARR.exists():
        raise FileNotFoundError(f"--append_anomalies requires existing store: {OUTPUT_ZARR}")
    print(f"\n=== APPEND MODE: adding anomaly arrays to {OUTPUT_ZARR} ===")

    # Need to re-ingest to get the per-profile T/Tclim arrays in memory.
    # (They are not stored in the zarr — only monthly means are.)
    source_dirs = _source_dirs(SOURCE_ROOT)
    bin_profiles   = defaultdict(lambda: defaultdict(list))
    bin_profiles_b = defaultdict(lambda: defaultdict(list))
    sources_seen   = set()

    all_work = []
    for sname, sdir in source_dirs:
        files = sorted(sdir.rglob('*.nc'))
        print(f"\n=== {sname} === {len(files)} files")
        all_work.extend((f, sname) for f in files)

    n_total = len(all_work)
    done    = 0
    print(f"\nIngesting {n_total} files with {_INGEST_WORKERS} workers (append mode)...")
    with concurrent.futures.ThreadPoolExecutor(max_workers=_INGEST_WORKERS) as pool:
        futures = {pool.submit(_ingest_file, item): item for item in all_work}
        for fut in concurrent.futures.as_completed(futures):
            source, fname, recs_fine, recs_coarse, err = fut.result()
            done += 1
            if err is not None:
                print(f"  [{done}/{n_total}] {fname}: skipping ({err})", flush=True)
                continue
            for bid, rec in recs_fine:
                bin_profiles[bid][source].append(rec)
            if recs_coarse:
                for bid_b, rec in recs_coarse:
                    bin_profiles_b[bid_b][source].append(rec)
            sources_seen.add(source)
            if done % 20 == 0 or done == n_total:
                print(f"  [{done}/{n_total}] ingested", flush=True)

    store = zarr.open_group(str(OUTPUT_ZARR), mode='a')
    _append_anomaly_arrays_from_store(store, 'fine',   bin_profiles)
    _append_anomaly_arrays_from_store(store, 'coarse', bin_profiles_b)
    zarr_size_mb = sum(f.stat().st_size for f in OUTPUT_ZARR.rglob('*') if f.is_file()) / 1e6
    print(f"\nAppend done — store size {zarr_size_mb:.1f} MB")
    import sys; sys.exit(0)

# ==============================================================================
# Full rebuild path
# ==============================================================================

# Reference depth grid
DEPTH_REF = None
for _, _d in _source_dirs(SOURCE_ROOT):
    for _f in sorted(_d.rglob('*.nc')):
        try: _ds = xr.open_dataset(_f)
        except Exception: continue
        if _ds.sizes.get('iDEPTH', 0) == N_DEPTH_REF:
            DEPTH_REF = _ds['depth'].values.astype(float).ravel()
            _ds.close(); break
        _ds.close()
    if DEPTH_REF is not None: break
if DEPTH_REF is None:
    raise RuntimeError(f"No post-NCEI file found with iDEPTH={N_DEPTH_REF}")
print(f"  Reference depth grid: {N_DEPTH_REF} levels, "
      f"{DEPTH_REF[0]:.0f}-{DEPTH_REF[-1]:.0f} m")

# Ingest
source_dirs = _source_dirs(SOURCE_ROOT)
bin_profiles   = defaultdict(lambda: defaultdict(list))
bin_profiles_b = defaultdict(lambda: defaultdict(list))
sources_seen   = set()

all_work = []
for sname, sdir in source_dirs:
    files = sorted(sdir.rglob('*.nc'))
    print(f"\n=== {sname} === {len(files)} files")
    all_work.extend((f, sname) for f in files)

n_total = len(all_work)
done    = 0
print(f"\nIngesting {n_total} files with {_INGEST_WORKERS} workers...")
with concurrent.futures.ThreadPoolExecutor(max_workers=_INGEST_WORKERS) as pool:
    futures = {pool.submit(_ingest_file, item): item for item in all_work}
    for fut in concurrent.futures.as_completed(futures):
        source, fname, recs_fine, recs_coarse, err = fut.result()
        done += 1
        if err is not None:
            print(f"  [{done}/{n_total}] {fname}: skipping ({err})", flush=True)
            continue
        for bid, rec in recs_fine:
            bin_profiles[bid][source].append(rec)
        if recs_coarse:
            for bid_b, rec in recs_coarse:
                bin_profiles_b[bid_b][source].append(rec)
        sources_seen.add(source)
        if done % 20 == 0 or done == n_total:
            print(f"  [{done}/{n_total}] ingested; "
                  f"{len(bin_profiles):,} fine / {len(bin_profiles_b):,} coarse bins so far",
                  flush=True)

sources_sorted = sorted(sources_seen)
print(f"\nSources ingested: {', '.join(sources_sorted)}")
print(f"  {len(bin_profiles):,} fine bins, {len(bin_profiles_b):,} coarse bins")

# Calendar span
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

# Aggregate
print(f"Aggregating {len(bin_profiles):,} fine bins...")
bins_out = aggregate_bins(bin_profiles, total_months)
print(f"Aggregating {len(bin_profiles_b):,} coarse bins...")
bins_out_b = aggregate_bins(bin_profiles_b, total_months)

# Coarse geodesic centers (full 2562 set)
_coarse_raw = np.genfromtxt(GEODESIC_FILE_B, delimiter=',')
_all_coarse_centers = {i + 1: (float(_coarse_raw[i, 0]), float(_coarse_raw[i, 1]))
                       for i in range(len(_coarse_raw))}
for k, v in _all_coarse_centers.items():
    bin_centers_b.setdefault(k, v)

def write_profiles_db(bins_fine, bins_coarse_dict):
    """Write per-profile records to SQLite for lazy lookup on bin click.

    Table: profiles(resolution TEXT, bin_id INT, src TEXT, year INT, month INT,
                    date INT, lon REAL, lat REAL, file TEXT, prof_idx INT, qual TEXT)
    Index on (resolution, bin_id) covers the only query pattern the app uses.
    """
    OUTPUT_DB = OUTPUT_CONTAINER / (OUTPUT_DIR.name + '_profiles.db')
    print(f"\nWriting profile records → {OUTPUT_DB} ...")
    if OUTPUT_DB.exists():
        OUTPUT_DB.unlink()
    con = sqlite3.connect(str(OUTPUT_DB))
    cur = con.cursor()
    cur.execute('''
        CREATE TABLE profiles (
            resolution TEXT NOT NULL,
            bin_id     INTEGER NOT NULL,
            src        TEXT NOT NULL,
            year       INTEGER NOT NULL,
            month      INTEGER NOT NULL,
            date       INTEGER NOT NULL,
            lon        REAL NOT NULL,
            lat        REAL NOT NULL,
            file       TEXT NOT NULL,
            prof_idx   INTEGER NOT NULL,
            qual       TEXT NOT NULL
        )
    ''')
    cur.execute('CREATE INDEX idx_bid ON profiles (resolution, bin_id)')

    def _rows(resolution, bins_dict):
        for bid, bd in bins_dict.items():
            for ym, m in bd['months'].items():
                for src, e in m['by_source'].items():
                    for s in e['sources']:
                        yield (resolution, bid, src, ym[0], ym[1],
                               s['date'], s['lon'], s['lat'],
                               s['file'], s['prof_idx'], s['qual'])

    con.executemany('INSERT INTO profiles VALUES (?,?,?,?,?,?,?,?,?,?,?)',
                    _rows('fine', bins_fine))
    con.executemany('INSERT INTO profiles VALUES (?,?,?,?,?,?,?,?,?,?,?)',
                    _rows('coarse', bins_coarse_dict))
    con.commit()
    con.close()
    db_size_mb = OUTPUT_DB.stat().st_size / 1e6
    print(f"  SQLite done — {db_size_mb:.1f} MB")


OUTPUT_CONTAINER.mkdir(parents=True, exist_ok=True)
write_zarr_and_json(
    bins_out, bins_out_b,
    bin_centers, bin_centers_b,
    DEPTH_REF.astype(np.float32),
    sources_sorted,
    (y0, m0) if all_ym else None,
    (y1, m1) if all_ym else None,
    total_months,
)
write_profiles_db(bins_out, bins_out_b)

print(f"\nDone. {len(bins_out):,} fine bins + {len(bins_out_b):,} coarse bins, "
      f"{len(sources_sorted)} sources.")
