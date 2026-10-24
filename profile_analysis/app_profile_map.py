"""
Profile map app — one flat map of the NCEI-processed profiles, three projections.

Replaces app_allsource_globe.py, app_allsource_anomaly_globe.py and
app_profile_scatter.py (the globe versions are kept at git tag globe-apps-2026-10-19).

Map
  Projection : World (natural earth), Arctic, Southern Ocean (polar stereographic,
               opening at map_common.POLAR_LAT_LIMIT)
  Show       : Bins     — geodesic bins colored by visit probability or by the
                          long-run T / S anomaly mean / std at a chosen depth
               Profiles — individual profiles as dots colored by source, within a
                          year window, over faint bin outlines
  Click a bin (or a dot) to load that bin's depth–time T / S panels below.

Reads the store built by build_allsource_bin_timeseries_with_anomalies.py
(zarr + JSON metadata + SQLite profile records); per-bin arrays and profile rows are
read when needed.

Run with:
    python app_profile_map.py [--data <store_dir>]
Then open http://127.0.0.1:8050
"""

import argparse
import concurrent.futures
import functools
import json
import math
import sqlite3
import threading
from pathlib import Path

import numpy as np
import plotly.colors as pcolors
import plotly.graph_objects as go
import zarr
from dash import Dash, dcc, html, Input, Output, State, dash_table, Patch, no_update, ctx

import map_common as mc
import subpolygons as sp

# ==============================================================================
# Config
# ==============================================================================

_OUTPUT_BASE = Path('/Users/brucel/ecco/yip/profile_data/z_profile_file_analysis/'
                    'allsource_bin_timeseries_with_anomalies')


def _latest_store(base):
    """Return the most recent timestamped store directory, falling back to base itself."""
    candidates = sorted(p for p in base.parent.glob(base.name + '_????????_??????') if p.is_dir())
    return candidates[-1] if candidates else base


_app_parser = argparse.ArgumentParser(description='Profile map app.')
_app_parser.add_argument('--data', default=None,
                         help='Store directory to load (default: the most recent '
                              'allsource_bin_timeseries_with_anomalies_* directory).')
_app_args, _unknown = _app_parser.parse_known_args()

_DATA_DIR     = Path(_app_args.data) if _app_args.data else _latest_store(_OUTPUT_BASE)
ZARR_STORE    = _DATA_DIR / (_DATA_DIR.name + '.zarr')
METADATA_JSON = _DATA_DIR / (_DATA_DIR.name + '_metadata.json')
PROFILES_DB   = _DATA_DIR / (_DATA_DIR.name + '_profiles.db')
_missing = [p.name for p in (ZARR_STORE, METADATA_JSON, PROFILES_DB) if not p.exists()]
if _missing:
    raise SystemExit(f'No complete data store at {_DATA_DIR} (missing {", ".join(_missing)}).\n'
                     'Pass --data <store folder>; build one with build_allsource_bin_timeseries_with_anomalies.py '
                     '(see README_maps.md, "What you need").')
print(f"Using data store: {_DATA_DIR}", flush=True)

SOURCE_ALL = 'ALL'
RES_FINE, RES_COARSE = 'fine', 'coarse'
SHOW_BINS, SHOW_PROFILES = 'bins', 'profiles'
TARGET_N = 20_000                 # Profiles view: at most this many dots on screen

# Bins view color modes
CMODE_PROB    = 'prob'
CMODE_T_MEAN  = 't_anom_mean'
CMODE_T_STD   = 't_anom_std'
CMODE_S_MEAN  = 's_anom_mean'
CMODE_S_STD   = 's_anom_std'

CMODE_OPTIONS = [
    {'label': ' Visit probability',  'value': CMODE_PROB},
    {'label': ' T anomaly mean',     'value': CMODE_T_MEAN},
    {'label': ' T anomaly std',      'value': CMODE_T_STD},
    {'label': ' S anomaly mean',     'value': CMODE_S_MEAN},
    {'label': ' S anomaly std',      'value': CMODE_S_STD},
]

# (colorscale, label, zmin, zmax, diverging) per mode
_CMODE_CS = {
    CMODE_PROB:   ('RdYlGn',  'P(≥1)',   0.0,  1.0,  False),
    CMODE_T_MEAN: ('RdBu_r',  'T anom mean (°C)',  None, None, True),
    CMODE_T_STD:  ('cividis', 'T anom std (°C)',   0.0,  None, False),
    CMODE_S_MEAN: ('RdBu_r',  'S anom mean (PSU)', None, None, True),
    CMODE_S_STD:  ('cividis', 'S anom std (PSU)',  0.0,  None, False),
}

PROB_KEYS = {'>=1': 1, '>=2': 2, '>=3': 3}
BASIS_DAY0        = 'day0'
BASIS_SINCE_FIRST = 'since_first'

# ==============================================================================
# Load metadata
# ==============================================================================

print("Loading metadata JSON...", flush=True)
with open(METADATA_JSON, 'r') as _f:
    _meta = json.load(_f)

depth              = np.array(_meta['depth_list'], dtype=np.float32)
SOURCES            = _meta['sources']
DATASET_FIRST_YM   = tuple(_meta['dataset_first_ym']) if _meta.get('dataset_first_ym') else None
DATASET_LAST_YM    = tuple(_meta['dataset_last_ym'])  if _meta.get('dataset_last_ym')  else None
DATASET_MONTHS     = _meta.get('total_months', 1)
# key names: 'bin_centers'; stores written before 2026-10-19 used 'bin_centres'
bin_centers        = {int(k): tuple(v) for k, v in (_meta.get('bin_centers') or _meta['bin_centres']).items()}
bin_centers_coarse = {int(k): tuple(v) for k, v in
                      (_meta.get('bin_centers_coarse') or _meta['bin_centres_coarse']).items()}
N_DEPTH            = len(depth)


def _parse_bin_meta(raw):
    months = {}
    for ym_str, m in raw['months'].items():
        y, mo = int(ym_str[:4]), int(ym_str[5:])
        months[(y, mo)] = {'n': m['n'],
                           'by_source': {src: {'n': e['n']} for src, e in m['by_source'].items()}}
    by_source = {}
    for src, roll in raw['by_source'].items():
        by_source[src] = {
            'hit_counts': {int(k): v for k, v in roll['hit_counts'].items()},
            'first_ym':   tuple(roll['first_ym']),
            'last_ym':    tuple(roll['last_ym']),
            'n_profiles': roll['n_profiles'],
            'has_estim':  roll['has_estim'],
        }
    return {
        'bin_id':     raw['bin_id'],
        'hit_counts': {int(k): v for k, v in raw['hit_counts'].items()},
        'first_ym':   tuple(raw['first_ym']),
        'last_ym':    tuple(raw['last_ym']),
        'prob_1plus': raw['prob_1plus'],
        'months':     months,
        'by_source':  by_source,
    }


bins        = {int(k): _parse_bin_meta(v) for k, v in _meta['bins_meta'].items()}
bins_coarse = {int(k): _parse_bin_meta(v) for k, v in _meta['bins_coarse_meta'].items()}
del _meta
print(f"  {len(bins):,} fine bins, {len(bins_coarse):,} coarse bins, "
      f"{len(SOURCES)} sources: {', '.join(SOURCES)}")

print("Opening zarr store...", flush=True)
_ZARR = zarr.open_group(str(ZARR_STORE), mode='r')

# Read-only SQLite profile records (bin tables, Profiles view)
_DB_CON = sqlite3.connect(f'file:{PROFILES_DB}?mode=ro', uri=True, check_same_thread=False)
_DB_CON.row_factory = sqlite3.Row
_DB_LOCK = threading.Lock()       # one connection shared by Dash's threads

# ==============================================================================
# Geodesic cells (both grids), coastline, bin outlines
# ==============================================================================

print("Building geodesic cells...", flush=True)
CELL_GEOJSON, CELL_IDS, CELL_RINGS = mc.build_cells(bin_centers)
COARSE_GEOJSON, COARSE_IDS, COARSE_RINGS = mc.build_cells(bin_centers_coarse)
_COAST_TRACE = mc.coastline_trace()
_BIN_OUTLINES = {RES_FINE:   mc.outline_trace(CELL_RINGS, set(bins)),
                 RES_COARSE: mc.outline_trace(COARSE_RINGS, set(bins_coarse))}
print(f"  {len(CELL_IDS):,} fine + {len(COARSE_IDS):,} coarse cells")


def _grid(res):
    if res == RES_COARSE:
        return bins_coarse, bin_centers_coarse, COARSE_GEOJSON, COARSE_RINGS
    return bins, bin_centers, CELL_GEOJSON, CELL_RINGS

# ==============================================================================
# Zarr / SQLite accessors
# ==============================================================================

def _query_bin_profiles(resolution, bid):
    """Profile record dicts for a bin, grouped by {(year, month): {src: [recs]}}."""
    with _DB_LOCK:
        rows = _DB_CON.execute(
            'SELECT src, year, month, date, lon, lat, file, prof_idx, qual '
            'FROM profiles WHERE resolution=? AND bin_id=?', (resolution, bid)).fetchall()
    result = {}
    for row in rows:
        result.setdefault((row['year'], row['month']), {}).setdefault(row['src'], []).append({
            'date': row['date'], 'lon': row['lon'], 'lat': row['lat'], 'file': row['file'],
            'prof_idx': row['prof_idx'], 'qual': row['qual'], 'source': row['src'],
        })
    return result


def _zarr_src_monthly(resolution, bid, src):
    """(ym_index (n_months,) YYYYMM, {key: (n_months, N_DEPTH) float32}) for one source."""
    node     = _ZARR[f"{resolution}/{bid}/by_source/{src}"]
    ym_index = np.asarray(_ZARR[f"{resolution}/{bid}/ym_index"], dtype=np.int32)
    arrays   = {k: np.asarray(node[k], dtype=np.float32) for k in ('T', 'S', 'Tclim', 'Sclim')}
    return ym_index, arrays


def _zarr_arrays(resolution, bid, ym, src):
    """One month's per-source T/S/Tclim/Sclim, each (N_DEPTH,) float32."""
    ym_int           = ym[0] * 100 + ym[1]
    ym_index, arrays = _zarr_src_monthly(resolution, bid, src)
    idx = int(np.searchsorted(ym_index, ym_int))
    if idx >= len(ym_index) or ym_index[idx] != ym_int:
        nan = np.full(N_DEPTH, np.nan, dtype=np.float32)
        return {k: nan for k in ('T', 'S', 'Tclim', 'Sclim')}
    return {k: arrays[k][idx] for k in ('T', 'S', 'Tclim', 'Sclim')}


def _zarr_means(resolution, bid, src=None):
    path = f"{resolution}/{bid}" if src is None else f"{resolution}/{bid}/by_source/{src}"
    node = _ZARR[path]
    return (np.asarray(node['T_mean'], dtype=np.float32),
            np.asarray(node['S_mean'], dtype=np.float32))


def _zarr_anom_stats(resolution, bid):
    """(T_anom_mean, T_anom_std, S_anom_mean, S_anom_std) as float32 arrays."""
    node = _ZARR[f"{resolution}/{bid}"]
    return (np.asarray(node['T_anom_mean'], dtype=np.float32),
            np.asarray(node['T_anom_std'],  dtype=np.float32),
            np.asarray(node['S_anom_mean'], dtype=np.float32),
            np.asarray(node['S_anom_std'],  dtype=np.float32))


# Anomaly stats (T mean, T std, S mean, S std) and valid counts (T, S) of every bin,
# per depth, for one grid and one source (or all sources pooled), kept in memory.
# They come from the store's summary tables (one read each, < 1 s); stores built
# before those existed fall back to reading every bin (~30 s, all sources only).
_ANOM_TABLES = {}         # (res, src) -> ({bin_id: (4, N_DEPTH) float32}, {bin_id: (2, N_DEPTH) int32})
_ANOM_LOCKS = {}
_NO_COUNTS = np.zeros((2, N_DEPTH), dtype=np.int32)


def _anom_tables(res, src=SOURCE_ALL):
    key = (res, src)
    with _ANOM_LOCKS.setdefault(key, threading.Lock()):
        if key not in _ANOM_TABLES:
            try:
                g = _ZARR[f"summary/{res}/{src}"]
                bids = np.asarray(g['bin_ids'])
                stats, n = np.asarray(g['stats'], dtype=np.float32), np.asarray(g['n'])
                tables = ({int(b): stats[i] for i, b in enumerate(bids)},
                          {int(b): n[i] for i, b in enumerate(bids)})
            except KeyError:
                tables = _anom_tables_per_bin(res) if src == SOURCE_ALL else ({}, {})
            _ANOM_TABLES[key] = tables
        return _ANOM_TABLES[key]


def _anom_tables_per_bin(res):
    """Older stores (no summary tables): read every bin's pooled stats and counts."""
    def one(bid):
        try:
            node = _ZARR[f"{res}/{bid}"]
            stats = np.stack(_zarr_anom_stats(res, bid))
        except Exception:
            return bid, None, _NO_COUNTS
        try:
            n = np.stack([np.asarray(node['T_anom_n']), np.asarray(node['S_anom_n'])])
        except Exception:
            n = _NO_COUNTS
        return bid, stats, n
    with concurrent.futures.ThreadPoolExecutor(max_workers=16) as ex:
        out = list(ex.map(one, sorted(_grid(res)[0])))
    return {b: st for b, st, _ in out if st is not None}, {b: n for b, _, n in out}


def _anom_rows(res, src=SOURCE_ALL):
    return _anom_tables(res, src)[0]


def _anom_counts(res, src=SOURCE_ALL):
    return _anom_tables(res, src)[1]


for _r in (RES_FINE, RES_COARSE):
    threading.Thread(target=_anom_tables, args=(_r,), daemon=True, name=f'anom-preload-{_r}').start()


def _preload_subpoly_ranges():
    """Read the single-profile anomalies once (all sources, T and S, fine grid) in the
    background (~1 min), so the first zoom to sub-polygons doesn't wait for it."""
    for var in (0, 1):
        _single_profile_values(RES_FINE, SOURCE_ALL, var)


@functools.lru_cache(maxsize=None)
def _cs(name):
    """A named color scale as an explicit list. Plotly's browser side does not know
    Python's reversed names (e.g. 'RdBu_r'), so a Patch must send the colors."""
    return [list(c) for c in pcolors.get_colorscale(name)]

# ==============================================================================
# Bin statistics
# ==============================================================================

def _months_between(ym0, ym1):
    return (ym1[0] - ym0[0]) * 12 + (ym1[1] - ym0[1]) + 1


def _pool_level_weighted(entries, key):
    num = np.zeros(N_DEPTH)
    den = np.zeros(N_DEPTH)
    for e in entries:
        v = np.asarray(e[key], dtype=float)
        mask = ~np.isnan(v)
        num[mask] += v[mask] * e['n']
        den[mask] += e['n']
    out = np.full(N_DEPTH, np.nan)
    nz = den > 0
    out[nz] = num[nz] / den[nz]
    return out


def bin_has_source(bin_data, source_sel):
    return source_sel == SOURCE_ALL or source_sel in bin_data['by_source']


def bin_probability(bin_data, source_sel, n_threshold, basis):
    if source_sel == SOURCE_ALL:
        hits, first_ym = bin_data['hit_counts'].get(n_threshold, 0), bin_data['first_ym']
    else:
        roll = bin_data['by_source'].get(source_sel)
        if roll is None:
            return 0.0
        hits, first_ym = roll['hit_counts'].get(n_threshold, 0), roll['first_ym']
    if basis == BASIS_SINCE_FIRST and first_ym and DATASET_LAST_YM:
        denom = _months_between(first_ym, DATASET_LAST_YM)
    else:
        denom = DATASET_MONTHS
    return hits / denom if denom > 0 else 0.0


def _visible_ids(source_sel, res=RES_FINE):
    b_dict, c_dict, _, _ = _grid(res)
    return [b for b in sorted(b_dict) if b in c_dict and bin_has_source(b_dict[b], source_sel)]


def resolve_bin(bin_data, source_sel, res=RES_FINE):
    bid = bin_data['bin_id']
    db_recs = _query_bin_profiles(res, bid)

    months = {}
    for ym, m in bin_data['months'].items():
        if source_sel != SOURCE_ALL and source_sel not in m['by_source']:
            continue
        src_list = list(m['by_source']) if source_sel == SOURCE_ALL else [source_sel]
        entries = []
        for src in src_list:
            if src not in m['by_source']:
                continue
            entry = dict(_zarr_arrays(res, bid, ym, src))
            entry['n'] = m['by_source'][src]['n']
            entries.append(entry)
        if not entries:
            continue
        recs = []
        for src in src_list:
            recs.extend(db_recs.get(ym, {}).get(src, []))
        months[ym] = {
            'T':       _pool_level_weighted(entries, 'T'),
            'S':       _pool_level_weighted(entries, 'S'),
            'Tclim':   _pool_level_weighted(entries, 'Tclim'),
            'Sclim':   _pool_level_weighted(entries, 'Sclim'),
            'sources': recs,
            'n':       sum(e['n'] for e in entries),
        }

    if source_sel == SOURCE_ALL:
        T_mean, S_mean = _zarr_means(res, bid)
        n_profiles = sum(r['n_profiles'] for r in bin_data['by_source'].values())
    else:
        T_mean, S_mean = _zarr_means(res, bid, source_sel)
        n_profiles = bin_data['by_source'][source_sel]['n_profiles']
    return {'bin_id': bid, 'months': months,
            'T_mean': np.asarray(T_mean, dtype=float), 'S_mean': np.asarray(S_mean, dtype=float),
            'n_profiles': n_profiles}

# ==============================================================================
# Bins view: color arrays and traces
#   trace order: [0] bins choropleth, [1] selected-bin ring, [2] coastline
# ==============================================================================

# ---- Color ranges --------------------------------------------------------------
# Fixed per selection (grid, source, variable, mean or std): computed over all bins and
# all depths, each (bin, depth) weighted by its number of profiles, so a color means the
# same at every depth and zoom. "Rescale to this depth" switches to the same statistic
# at the chosen depth only (more detail where anomalies are small, e.g. at depth); it
# turns itself off at RESCALE_RESET_DEPTH and above, and with "Reset view".
RANGE_PCTL = 95.0                # mean modes: ±P95 of |value|; std modes: 0 .. P95
RESCALE_RESET_DEPTH = 200.0      # m: "rescale to this depth" switches off at or above this


def _wpct(values, weights, pct):
    """Weighted percentile (values, weights 1-D; NaNs ignored)."""
    ok = np.isfinite(values) & (weights > 0)
    if not ok.any():
        return np.nan
    v, w = values[ok], weights[ok]
    o = np.argsort(v); v, w = v[o], w[o]
    c = np.cumsum(w)
    return float(np.interp(pct / 100.0 * c[-1], c, v))


@functools.lru_cache(maxsize=512)
def _color_range(res, src, k, depth_idx=None):
    """Upper end of the bins' color range for stat row k (0 T mean, 1 T std, 2 S mean,
    3 S std): P95 of |value| over (bins x depths), weighted by profile counts; one depth
    only if depth_idx is given."""
    rows, counts = _anom_tables(res, src)
    ids = [b for b in rows if b in counts]
    if not ids:
        return 1.0
    S = np.stack([rows[b] for b in ids])                      # (nbins, 4, ndepth)
    N = np.stack([counts[b] for b in ids]).astype(float)       # (nbins, 2, ndepth)
    var = 0 if k in (0, 1) else 1
    x, w = S[:, k, :], N[:, var, :]
    if depth_idx is not None:
        x, w = x[:, depth_idx], w[:, depth_idx]
    hi = _wpct(np.abs(x.ravel()), w.ravel(), RANGE_PCTL)
    return hi if np.isfinite(hi) and hi > 1e-9 else 1.0


def _subpoly_range(res, src, var, depth_idx=None):
    """Upper end of the sub-polygons' color range: P95 of |single-profile anomaly|
    (var 0 T, 1 S); falls back to the bins' range where there are no values."""
    single = _single_profile_range(res, src, var, depth_idx)
    return single if np.isfinite(single) and single > 1e-9 else _color_range(res, src, 2 * var, depth_idx)


_SP_LOCK = threading.Lock()


def _single_profile_values(res, src, var):
    with _SP_LOCK:                                   # one reader at a time; repeats come from the cache
        return _single_profile_values_cached(res, src, var)


@functools.lru_cache(maxsize=256)
def _single_profile_values_cached(res, src, var):
    """|single-profile anomaly| values and weights, per depth, from every bin's stored
    quantiles (var 0 T, 1 S): {depth_idx: (values, weights)}. Read once per selection
    (all bins, in parallel), then kept."""
    ids = sorted(_anom_tables(res, src)[0])
    def one(bid):
        q = _anom_quantiles(res, bid, src)
        return None if q is None else (q[2 * var], q[2 * var + 1])
    with concurrent.futures.ThreadPoolExecutor(max_workers=16) as ex:
        got = [g for g in ex.map(one, ids) if g is not None]
    out = {}
    for d in range(N_DEPTH):
        v, w = [], []
        for q, n in got:
            kk = min(int(n[d]), q.shape[0])
            if kk:
                v.append(np.abs(q[:kk, d])); w.append(np.full(kk, n[d] / kk))
        out[d] = (np.concatenate(v), np.concatenate(w)) if v else (np.empty(0), np.empty(0))
    return out


def _single_profile_range(res, src, var, depth_idx=None):
    vals = _single_profile_values(res, src, var)
    if depth_idx is not None:
        v, w = vals[depth_idx]
    else:
        v = np.concatenate([a for a, _ in vals.values()]); w = np.concatenate([b for _, b in vals.values()])
    return _wpct(v, w, RANGE_PCTL)


def _anom_arrays_for_display(ids_clean, cmode, depth_idx, res, src=SOURCE_ALL, rescale=False):
    """(z_values, colorscale, zmin, zmax, label, depth_m) for the anomaly modes."""
    colorscale, label, zmin_fixed, zmax_fixed, diverging = _CMODE_CS[cmode]
    k = {CMODE_T_MEAN: 0, CMODE_T_STD: 1, CMODE_S_MEAN: 2}.get(cmode, 3)   # row in _anom_rows
    rows = _anom_rows(res, src)
    # None (no fill) where the source has no data at this depth: drawing those bins
    # white would read as a zero anomaly
    z_vals = []
    for bid in ids_clean:
        r = rows.get(bid)
        z_vals.append(float(r[k, depth_idx]) if r is not None and not np.isnan(r[k, depth_idx]) else None)
    hi = _color_range(res, src, k, depth_idx if rescale else None)
    zmin, zmax = (-hi, hi) if diverging else (0.0, hi)
    return z_vals, _cs(colorscale), zmin, zmax, label, float(depth[depth_idx])


def _bin_arrays(source_sel, n_threshold, basis, cmode, depth_idx, res=RES_FINE, rescale=False):
    b_dict, c_dict, _, _ = _grid(res)
    ids_clean  = _visible_ids(source_sel, res)
    prob_label = next(k for k, v in PROB_KEYS.items() if v == n_threshold)
    src_note   = 'All sources' if source_sel == SOURCE_ALL else source_sel
    probs = [bin_probability(b_dict[b], source_sel, n_threshold, basis) for b in ids_clean]

    if cmode == CMODE_PROB:
        hover = [f'Bin {b} — {src_note}<br>lon={c_dict[b][0]:.2f}, lat={c_dict[b][1]:.2f}'
                 f'<br>P({prob_label} prof/month) = {p:.3f}' for b, p in zip(ids_clean, probs)]
        return dict(ids=ids_clean, z=probs, hover=hover,
                    colorscale=_cs('RdYlGn'), zmin=0.0, zmax=1.0,
                    cb_title=f'P({prob_label}<br>prof/month)')        # same wording as the hover

    z_vals, colorscale, zmin, zmax, label, depth_m = _anom_arrays_for_display(
        ids_clean, cmode, depth_idx, res, source_sel, rescale)
    counts = _anom_counts(res, source_sel)
    kv = 0 if cmode in (CMODE_T_MEAN, CMODE_T_STD) else 1
    hover = [f'Bin {b} — {src_note}<br>lon={c_dict[b][0]:.2f}, lat={c_dict[b][1]:.2f}'
             + (f'<br>{label} @ {depth_m:.0f} m = {z:.3f}' if z is not None else f'<br>no data at {depth_m:.0f} m')
             + f'<br>profiles with data at this depth: {int(counts.get(b, _NO_COUNTS)[kv, depth_idx]):,}'
             + f'<br>P({prob_label}) = {p:.3f}'
             for b, z, p in zip(ids_clean, z_vals, probs)]
    return dict(ids=ids_clean, z=z_vals, hover=hover, colorscale=colorscale, zmin=zmin, zmax=zmax,
                cb_title=label.replace(' (', '<br>('))         # two lines: fits beside the edge bar


def _ring_lonlat(bid, res):
    ring = _grid(res)[3].get(bid) if bid is not None else None
    if ring is None:
        return [], []
    return ring[:, 0].tolist(), ring[:, 1].tolist()


def _selected_trace(selected_bin_id, res=RES_FINE):
    lon, lat = _ring_lonlat(selected_bin_id, res)
    return go.Scattergeo(lon=lon, lat=lat, mode='lines',
                         line=dict(color='rgba(255,255,255,0.95)', width=2.5),
                         hoverinfo='skip', showlegend=False, name='selection')


def _make_choropleth(a, geojson):
    return go.Choropleth(
        geojson=geojson, locations=a['ids'], z=a['z'], featureidkey='id',
        colorscale=a['colorscale'], zmin=a['zmin'], zmax=a['zmax'],
        marker=dict(line=dict(color='rgba(40,40,40,0.5)', width=0.3)),
        colorbar=dict(title=dict(text=a['cb_title'], side='top', font=dict(size=11, color='#ddd')),
                      tickfont=dict(size=11, color='#ccc'), **CBAR_FILL),
        text=a['hover'], hoverinfo='text', customdata=a['ids'],
    )


# ---- Inset outlines (anomaly modes): each bin's ring shrunk EDGE_INSET toward its
# center, so no two bins share an edge.
#   color: the bin's anomaly std at the chosen depth (mean modes; gray in std modes,
#           where the fill already is the std)
#   width : zoomed out, the number of profiles with data at that depth (linear, capped at
#           the busiest 2% of bins; none there = no outline), so dense areas stand out;
#           while sub-polygons show, one thin uniform width (only the color matters)
EDGE_INSET = 0.08
EDGE_WIDTH_RANGE = (0.2, 3.0)    # px, fewest → most profiles, at the opening zoom (×√2 per zoom doubling)
EDGE_WIDTH_SUBPOLY = 1.5         # px, uniform, while sub-polygons show
EDGE_WIDTH_EMPTY = 0.4           # px: bins with the source's profiles but none at this depth
EDGE_COLOR_EMPTY = 'rgba(160,160,160,0.7)'
_STD_CS = 'cividis'

# Bins view layout: the map on the left, the color bars in a column on the right.
# Bins view layout: the map on the left; the color bars stacked in one column on the
# right, fill on top and edge below (side by side, their titles overlapped at some widths).
# Bins view layout: the edge (std) color bar on the left, the fill color bar on the
# right, the map between them; both bars tall, titles on top.
MAP_DOMAIN_X = [0.075, 0.925]
CBAR_FILL = dict(x=0.935, xanchor='left', y=0.47, yanchor='middle', len=0.8, thickness=10)
CBAR_EDGE = dict(x=0.005, xanchor='left', y=0.47, yanchor='middle', len=0.8, thickness=10)


def _inset_geojson(rings, centers):
    feats = [{'type': 'Feature', 'id': int(b),
              'geometry': {'type': 'Polygon',
                           'coordinates': [np.round(sp.inset_ring(rings[b], centers[b], EDGE_INSET), 3).tolist()]},
              'properties': {}} for b in sorted(rings)]
    return {'type': 'FeatureCollection', 'features': feats}


print("Building inset bin outlines...", flush=True)
_INSET_GEOJSON = {RES_FINE:   _inset_geojson(CELL_RINGS, bin_centers),
                  RES_COARSE: _inset_geojson(COARSE_RINGS, bin_centers_coarse)}


def _edges(ids, cmode, depth_idx, res, bucket=1, uniform=False, src=SOURCE_ALL, rescale=False):
    """Outline colors / widths for the bins in ids, or None when not drawn (probability)."""
    if cmode == CMODE_PROB:
        return None
    t_var = cmode in (CMODE_T_MEAN, CMODE_T_STD)
    counts = _anom_counts(res, src)
    n = np.array([int(counts.get(b, _NO_COUNTS)[0 if t_var else 1, depth_idx]) for b in ids])
    pos = n[n > 0]
    nmax = max(float(np.percentile(pos, 98)) if pos.size else 1.0, 1.0)
    # bins where the source has profiles, but none at this depth: a thin gray outline,
    # so the source's footprint stays visible (their fill is left empty)
    empty_w = EDGE_WIDTH_EMPTY * (1.0 if uniform else np.sqrt(bucket))
    if uniform:
        widths = np.where(n > 0, EDGE_WIDTH_SUBPOLY, empty_w)
    else:
        lo, hi = EDGE_WIDTH_RANGE
        widths = np.where(n > 0, (lo + (hi - lo) * np.clip(n / nmax, 0, 1)) * np.sqrt(bucket), empty_w)
    out = dict(widths=np.round(widths, 2).tolist(), nmax=int(nmax), cmax=None, title='',
               uniform=uniform, n_with=int((n > 0).sum()), n_bins=len(ids),
               src=src, depth_m=float(depth[depth_idx]))
    if cmode in (CMODE_T_MEAN, CMODE_S_MEAN):
        rows = _anom_rows(res, src)
        k = 1 if t_var else 3
        std = np.array([float(rows[b][k, depth_idx]) if rows.get(b) is not None else np.nan for b in ids])
        cmax = _color_range(res, src, k, depth_idx if rescale else None)
        frac = np.clip(np.nan_to_num(std, nan=0.0) / cmax, 0, 1)
        colors = pcolors.sample_colorscale(_STD_CS, frac.tolist())
        out.update(cmax=cmax, title=('edge color<br>' + ('T' if t_var else 'S') + ' anom std<br>'
                                     + ('(°C)' if t_var else '(PSU)')))
    else:
        colors = ['rgba(225,225,225,0.85)'] * len(ids)
    out['colors'] = [c if k > 0 else EDGE_COLOR_EMPTY for c, k in zip(colors, n)]
    return out


def _edge_note(e):
    if not e:
        return ''
    src = 'All sources' if e['src'] == SOURCE_ALL else e['src']
    head = (f'{src}: data at {e["depth_m"]:.0f} m in {e["n_with"]:,} of {e["n_bins"]:,} bins'
            + (' (gray outline = profiles, but none at this depth)' if e['n_with'] < e['n_bins'] else ''))
    color_note = '; color = anomaly std' if e['cmax'] is not None else ''
    if e['uniform']:
        return head + '<br>Edges: uniform width while sub-polygons show' + color_note
    return head + (f'<br>Edges: width ∝ profiles with data at this depth (thickest = {e["nmax"]:,} or more)'
                   + color_note)


def _edge_traces(ids, cmode, depth_idx, res, show_edges, bucket=1, src=SOURCE_ALL, rescale=False):
    """[inset outline choropleth, edge color-bar carrier]; hidden when not drawn."""
    e = _edges(ids, cmode, depth_idx, res, bucket, src=src, rescale=rescale) if show_edges else None
    outline = go.Choropleth(
        geojson=_INSET_GEOJSON[res], locations=ids, z=[0] * len(ids), featureidkey='id',
        colorscale=[[0, 'rgba(0,0,0,0)'], [1, 'rgba(0,0,0,0)']], showscale=False,
        marker=dict(line=dict(color=e['colors'] if e else 'rgba(0,0,0,0)',
                              width=e['widths'] if e else 0)),
        hoverinfo='skip', visible=e is not None, name='_edges')
    cbar = go.Scattergeo(
        lon=[None], lat=[None], mode='markers', hoverinfo='skip', showlegend=False, name='_edge_cbar',
        visible=bool(e and e['cmax'] is not None),
        marker=dict(color=[0.0], colorscale=_cs(_STD_CS), cmin=0.0,
                    cmax=(e['cmax'] if e and e['cmax'] is not None else 1.0), showscale=True,
                    colorbar=dict(title=dict(text=e['title'] if e else '', side='top',
                                             font=dict(size=11, color='#ddd')),
                                  tickfont=dict(size=11, color='#ccc'), **CBAR_EDGE)))
    return [outline, cbar]


# ---- Sub-polygons (zoomed in, anomaly-mean modes): each bin in view split into
# min(n, K) cells, one per profile when n <= K, else one per (i+0.5)/K quantile of
# its anomalies at the chosen depth; values low to high from south to north.
SUBPOLY_VIEW_DEG = 15.0          # they appear once the view is less than this many degrees of latitude
                                 # tall: the same for every source and projection (World ≈ zoom 12,
                                 # polar ≈ zoom 5.6; ≤ ~150 fine bins in view even in the densest areas)
SUBPOLY_MAX_BINS = 150           # safety limit: never more bins than this at once


@functools.lru_cache(maxsize=4000)
def _anom_quantiles(res, bid, src=SOURCE_ALL):
    """(T_q, T_n, S_q, S_n) of one source (or all pooled) from the store, or None if
    it has no such quantile arrays."""
    try:
        node = _ZARR[f"{res}/{bid}" if src == SOURCE_ALL else f"{res}/{bid}/by_source/{src}"]
        return (np.asarray(node['T_anom_q'], dtype=np.float32), np.asarray(node['T_anom_n']),
                np.asarray(node['S_anom_q'], dtype=np.float32), np.asarray(node['S_anom_n']))
    except Exception:
        return None


@functools.lru_cache(maxsize=4000)
def _sub_cells(res, bid, k):
    """Sub-polygon rings of a bin, in south-to-north order (cached)."""
    ring = _grid(res)[3][bid]
    cells = sp.subdivide(ring, _grid(res)[1][bid], k, seed=bid)
    return tuple(cells[i] for i in sp.order_south_to_north(cells))


def _in_view(lon, lat, state, margin=1.15):
    """Boolean mask: which lon / lat points are inside the current view (with a margin)."""
    lon, lat = np.asarray(lon, dtype=float), np.asarray(lat, dtype=float)
    z = mc.zoom_factor(state)
    if state['proj'] == mc.PROJ_WORLD:
        dlon = (lon - state['lon'] + 180.0) % 360.0 - 180.0
        return (np.abs(dlon) <= margin * 180.0 / z) & (np.abs(lat - state['lat']) <= margin * 90.0 / z)
    # stereographic: radius ∝ tan(colatitude / 2); the opening view spans the
    # POLAR_LAT_LIMIT circle, so at zoom z it spans colatitude 2·atan(tan(c0/2)/z)
    c0 = np.radians(90.0 - mc.POLAR_LAT_LIMIT)
    radius = margin * 2.0 * np.arctan(np.tan(c0 / 2.0) / z) * 1.42       # square map: corners
    p0, l0 = np.radians(state['lat']), np.radians(state['lon'])
    p, l = np.radians(lat), np.radians(lon)
    cosd = np.sin(p0) * np.sin(p) + np.cos(p0) * np.cos(p) * np.cos(l - l0)
    return np.arccos(np.clip(cosd, -1.0, 1.0)) <= radius


def _view_height_deg(state):
    """Degrees of latitude the current view spans, top to bottom."""
    z = mc.zoom_factor(state)
    if state['proj'] == mc.PROJ_WORLD:
        return 180.0 / z
    c0 = np.radians(90.0 - mc.POLAR_LAT_LIMIT)
    return 2.0 * np.degrees(2.0 * np.arctan(np.tan(c0 / 2.0) / z))


def _bins_in_view(ids, state, res):
    c_dict = _grid(res)[1]
    keep = _in_view([c_dict[b][0] for b in ids], [c_dict[b][1] for b in ids], state)
    return [b for b, k in zip(ids, keep) if k]


def _subpoly_layer(state, source_sel, cmode, depth_idx, res, rescale=False):
    """(geojson, ids, z, hover, zmin, zmax, colorscale, note) for the sub-polygon trace."""
    empty = ({'type': 'FeatureCollection', 'features': []}, [], [], [], 0.0, 1.0, _cs('RdBu_r'), '')
    if cmode not in (CMODE_T_MEAN, CMODE_S_MEAN) or not state or _view_height_deg(state) > SUBPOLY_VIEW_DEG:
        return empty
    ids = _bins_in_view(_visible_ids(source_sel, res), state, res)
    if len(ids) > SUBPOLY_MAX_BINS:
        return empty[:-1] + (f'{len(ids)} bins in view (limit {SUBPOLY_MAX_BINS}) — zoom in a little more '
                             f'for sub-polygons',)
    letter = 'T' if cmode == CMODE_T_MEAN else 'S'
    units = '°C' if letter == 'T' else 'PSU'
    depth_m = float(depth[depth_idx])
    src_note = 'all sources' if source_sel == SOURCE_ALL else source_sel
    feats, locs, z, hover = [], [], [], []
    var = 0 if letter == 'T' else 2
    for bid in ids:
        qa = _anom_quantiles(res, bid, source_sel)
        if qa is None:
            return empty[:-1] + ('this store has no ' + ('' if source_sel == SOURCE_ALL else 'per-source ')
                                 + 'quantile arrays (update it with --append_anomalies)',)
        q, n = qa[var], int(qa[var + 1][depth_idx])
        k = min(n, q.shape[0])
        if k == 0:
            continue
        cells = _sub_cells(res, bid, k)
        vals = q[:k, depth_idx]
        pick = np.round(np.linspace(0, k - 1, len(cells))).astype(int)   # if a cell failed
        for j, (cell, i) in enumerate(zip(cells, pick)):
            fid = f'{bid}:{j}'
            feats.append({'type': 'Feature', 'id': fid, 'properties': {},
                          'geometry': {'type': 'Polygon', 'coordinates': [cell.tolist()]}})
            locs.append(fid); z.append(float(vals[i]))
            if n <= q.shape[0]:
                what = f'Profile {i + 1} of {n}'
            else:
                what = f'Percentile {100 * (i + 0.5) / k:.1f} of {n:,} profiles'
            hover.append(f'Bin {bid} · {depth_m:.0f} m · {src_note}<br>{what}<br>'
                         f'{letter} − climatology = {vals[i]:+.3f} {units}')
    # Single profiles vary far more than bin means, so the cells get their own range,
    # fixed like the bins' (all bins and depths, or this depth when rescaled).
    zabs = _subpoly_range(res, source_sel, 0 if letter == 'T' else 1, depth_idx if rescale else None)
    return ({'type': 'FeatureCollection', 'features': feats}, locs, z, hover, -zabs, zabs, _cs('RdBu_r'),
            f'sub-polygons: {len(ids)} bins' if locs else '')


def _subpoly_cbar_title(cmode):
    letter = 'T' if cmode == CMODE_T_MEAN else 'S'
    return f'{letter} − clim<br>single profiles<br>' + ('(°C)' if letter == 'T' else '(PSU)')


def _subpoly_trace(layer=None):
    gj, locs, z, hover, zmin, zmax, cs, _ = layer or _subpoly_layer(None, SOURCE_ALL, CMODE_PROB, 0, RES_FINE)
    return go.Choropleth(geojson=gj, locations=locs, z=z, featureidkey='id', colorscale=cs,
                         zmin=zmin, zmax=zmax, showscale=bool(locs), text=hover, hoverinfo='text',
                         colorbar=dict(title=dict(text='', side='top', font=dict(size=11, color='#ddd')),
                                       tickfont=dict(size=11, color='#ccc'), **CBAR_FILL),
                         customdata=[int(f.split(':')[0]) for f in locs],
                         marker=dict(line=dict(color='rgba(20,20,20,0.6)', width=0.3)), name='_subpoly')


def bins_traces(source_sel, selected_bin_id, n_threshold, basis, res, cmode, depth_idx,
                show_edges=True, state=None, rescale=False):
    """[0] bins, [1] selected-bin ring, [2] coastline, [3] sub-polygons,
    [4] inset outlines, [5] their color bar."""
    a = _bin_arrays(source_sel, n_threshold, basis, cmode, depth_idx, res, rescale)
    return ([_make_choropleth(a, _grid(res)[2]), _selected_trace(selected_bin_id, res), _COAST_TRACE,
             _subpoly_trace()]
            + _edge_traces(a['ids'], cmode, depth_idx, res, show_edges, mc.zoom_bucket(state), src=source_sel,
                           rescale=rescale))

threading.Thread(target=_preload_subpoly_ranges, daemon=True, name='subpoly-range-preload').start()

# ==============================================================================
# Profiles view: dots
#   trace order: [0] bin outlines, [1..n] dots (n = _n_dot_traces), [n+1] star, [n+2] coastline
# ==============================================================================

_PALETTE = ['#636EFA', '#EF553B', '#00CC96', '#AB63FA', '#FFA15A',
            '#19D3F3', '#FF6692', '#B6E880', '#FF97FF', '#FECB52']
_SRC_COLOR = {src: _PALETTE[i % len(_PALETTE)] for i, src in enumerate(SOURCES)}


def _decimate(lats, lons, target_n=TARGET_N):
    """Boolean mask keeping ≤target_n spatially spread points.

    Cell size comes from the data's bounding box, so a tight regional query keeps
    more points than a global one.
    """
    n = len(lats)
    if n <= target_n:
        return np.ones(n, dtype=bool)
    lat_km   = float(lats.max() - lats.min()) * 111.0
    lon_km   = float(lons.max() - lons.min()) * 111.0 * max(0.1, math.cos(math.radians(float(np.mean(lats)))))
    cell_deg = max(0.01, math.sqrt(max(lat_km * lon_km, 1.0) / target_n) / 111.0)
    lb = np.floor(lats / cell_deg).astype(np.int64)
    bb = np.floor(lons / cell_deg).astype(np.int64)
    # Shuffle so the kept point per cell is random, not always the first in db order
    perm = np.random.default_rng(seed=0).permutation(n)
    _, first_in_perm = np.unique(lb[perm] * 1_000_000 + bb[perm], return_index=True)
    mask = np.zeros(n, dtype=bool)
    mask[perm[first_in_perm]] = True
    return mask


# Every profile's position, bin, source and year-month, per grid, read once from
# SQLite into compact arrays (~100 MB per grid); the Profiles view then filters in
# memory instead of re-reading ~6M rows on each source / year change. The fine grid
# is loaded in the background at start-up; the coarse grid when first needed.
_SRC_CODE = {s: i for i, s in enumerate(SOURCES)}
_DOTS = {}
_DOTS_LOCKS = {}


def _dots_table(res):
    with _DOTS_LOCKS.setdefault(res, threading.Lock()):
        if res not in _DOTS:
            con = sqlite3.connect(f'file:{PROFILES_DB}?mode=ro', uri=True)   # own connection, plain tuples
            cur = con.execute('SELECT lat, lon, bin_id, src, year * 100 + month FROM profiles '
                              'WHERE resolution=?', (res,))
            parts = {k: [] for k in ('lat', 'lon', 'bin', 'src', 'ym')}
            while True:
                rows = cur.fetchmany(500_000)
                if not rows:
                    break
                lat, lon, bid, src, ym = zip(*rows)
                parts['lat'].append(np.array(lat, dtype=np.float32))
                parts['lon'].append(np.array(lon, dtype=np.float32))
                parts['bin'].append(np.array(bid, dtype=np.int32))
                parts['src'].append(np.array([_SRC_CODE.get(s, -1) for s in src], dtype=np.int8))
                parts['ym'].append(np.array(ym, dtype=np.int32))
            con.close()
            _DOTS[res] = {k: (np.concatenate(v) if v else np.empty(0)) for k, v in parts.items()}
        return _DOTS[res]


threading.Thread(target=_dots_table, args=(RES_FINE,), daemon=True, name='dots-preload-fine').start()


def query_dots(source, res, year_min, year_max, view=None):
    """Profile locations in a year window (and, if view is given, in that map view),
    spatially decimated; None if there are none."""
    t = _dots_table(res)
    keep = (t['ym'] >= year_min * 100 + 1) & (t['ym'] <= year_max * 100 + 12)
    if source != SOURCE_ALL:
        keep &= t['src'] == _SRC_CODE.get(source, -2)
    if view is not None:
        idx = np.flatnonzero(keep)
        keep[idx[~_in_view(t['lon'][idx], t['lat'][idx], view)]] = False
    idx = np.flatnonzero(keep)
    if idx.size == 0:
        return None
    sub = _decimate(t['lat'][idx], t['lon'][idx])
    idx = idx[sub]
    src_names = np.array(SOURCES)
    return {'lats': t['lat'][idx], 'lons': t['lon'][idx], 'bin_ids': t['bin'][idx],
            'srcs': src_names[t['src'][idx]], 'dates': t['ym'][idx],
            'n_shown': int(idx.size), 'n_total': int(keep.sum())}


def _n_dot_traces(source_sel):
    return len(SOURCES) if source_sel == SOURCE_ALL else 1


def _dot_traces(data, source_sel):
    """Exactly _n_dot_traces(source_sel) traces (empty for sources with no profiles in
    the window), so the star that follows them is at a known index."""
    def _empty(name='_empty'):
        return go.Scattergeo(lat=[], lon=[], mode='markers', showlegend=False, name=name)

    def _dots(mask, src, color):
        d = {k: data[k][mask] for k in ('lats', 'lons', 'bin_ids', 'dates', 'srcs')}
        hover = [f'{s}<br>{dt // 100}-{dt % 100:02d}<br>lon={lo:.2f}, lat={la:.2f}'
                 for s, dt, lo, la in zip(d['srcs'], d['dates'], d['lons'], d['lats'])]
        return go.Scattergeo(lat=d['lats'].tolist(), lon=d['lons'].tolist(), mode='markers',
                             marker=dict(size=3, color=color, opacity=0.7),
                             name=src, text=hover, hoverinfo='text', customdata=d['bin_ids'].tolist())

    if data is None:
        return [_empty() for _ in range(_n_dot_traces(source_sel))]
    if source_sel == SOURCE_ALL:
        out = []
        for src in SOURCES:
            mask = data['srcs'] == src
            out.append(_dots(mask, src, _SRC_COLOR[src]) if mask.any() else _empty(src))
        return out
    return [_dots(np.ones(len(data['lats']), dtype=bool), source_sel,
                  _SRC_COLOR.get(source_sel, '#88ccff'))]


def _star_trace(sel_lonlat):
    if not sel_lonlat:
        return go.Scattergeo(lat=[], lon=[], mode='markers', hoverinfo='skip',
                             showlegend=False, name='_sel')
    return go.Scattergeo(lat=[sel_lonlat['lat']], lon=[sel_lonlat['lon']], mode='markers',
                         marker=dict(size=14, symbol='star', color='white', line=dict(color='black', width=1)),
                         hoverinfo='skip', showlegend=False, name='_sel')


_LEGEND_STYLE = dict(bgcolor='rgba(20,20,20,0.75)', font=dict(color='#ccc', size=11),
                     x=0.01, y=0.01, xanchor='left', yanchor='bottom')


def profiles_traces(data, source_sel, sel_lonlat, res):
    return ([_BIN_OUTLINES[res]] + _dot_traces(data, source_sel)
            + [_star_trace(sel_lonlat), _COAST_TRACE])

# ==============================================================================
# The map figure
# ==============================================================================

def make_map(proj, show, res, source_sel, uirevision, *, selected_bin=None, sel_lonlat=None,
             n_threshold=1, basis=BASIS_DAY0, cmode=CMODE_PROB, depth_idx=0, dots=None,
             show_edges=True, state=None, rescale=False):
    if show == SHOW_PROFILES:
        fig = go.Figure(profiles_traces(dots, source_sel, sel_lonlat, res))
        fig.update_layout(legend=_LEGEND_STYLE)
    else:
        fig = go.Figure(bins_traces(source_sel, selected_bin, n_threshold, basis, res, cmode, depth_idx,
                                    show_edges=show_edges, state=state, rescale=rescale))
        e = (_edges(_visible_ids(source_sel, res), cmode, depth_idx, res, src=source_sel, rescale=rescale)
             if show_edges else None)
        fig.add_annotation(text=_edge_note(e), x=MAP_DOMAIN_X[0], y=0.0, xref='paper', yref='paper', align='left',
                           xanchor='left', yanchor='bottom', showarrow=False,
                           font=dict(size=11, color='#ccc'), bgcolor='rgba(20,20,20,0.7)')
    mc.apply_projection(fig, proj, uirevision)
    if show != SHOW_PROFILES:
        fig.update_geos(domain=dict(x=MAP_DOMAIN_X, y=[0, 1]))
    return fig

# ==============================================================================
# Depth–time panels
# ==============================================================================

def make_zt_figure(bin_data, var, depth):
    months_dict = bin_data['months']
    sorted_yms  = sorted(months_dict)
    if not sorted_yms:
        return _empty_fig('No data for this source in this bin')
    x = [f'{ym[0]}-{ym[1]:02d}-15' for ym in sorted_yms]
    base_letter = var[0]
    units = '°C' if base_letter == 'T' else 'PSU'
    if var in ('T_minus_clim', 'S_minus_clim'):
        matrix = np.stack([months_dict[ym][base_letter] - months_dict[ym][base_letter + 'clim']
                           for ym in sorted_yms], axis=0).T
        title_var = f'{base_letter} − climatology'
    elif var in ('T_anom', 'S_anom'):
        matrix = np.stack([months_dict[ym][base_letter] for ym in sorted_yms], axis=0).T
        matrix = matrix - bin_data[f'{base_letter}_mean'][:, np.newaxis]
        title_var = f'{base_letter} − bin time-mean'
    else:
        matrix = None
    if matrix is not None:
        colorscale = 'RdBu_r'
        abs_max = np.nanpercentile(np.abs(matrix[np.isfinite(matrix)]), 97) if np.any(np.isfinite(matrix)) else 1
        zmin, zmax = -abs_max, abs_max
    else:
        matrix = np.stack([months_dict[ym][var] for ym in sorted_yms], axis=0).T
        valid  = matrix[np.isfinite(matrix)]
        colorscale = 'Thermal' if base_letter == 'T' else 'Haline'
        lo, hi = (-2, 30) if base_letter == 'T' else (30, 38)
        zmin = np.nanpercentile(valid, 2)  if len(valid) else lo
        zmax = np.nanpercentile(valid, 98) if len(valid) else hi
        title_var = base_letter + (' (clim)' if var.endswith('clim') else '')

    fig = go.Figure(go.Heatmap(
        x=x, y=depth, z=matrix, colorscale=colorscale, zmin=zmin, zmax=zmax,
        colorbar=dict(title=units, thickness=12),
        hovertemplate='%{x|%b %Y}<br>Depth: %{y} m<br>Value: %{z:.3f}<extra></extra>',
    ))
    bid = bin_data['bin_id']
    lon, lat = (bin_centers.get(bid) or bin_centers_coarse.get(bid) or (np.nan, np.nan))
    fig.update_layout(
        title=dict(text=f'Bin {bid} ({lon:.1f}°, {lat:.1f}°) — {title_var} (monthly bin average)',
                   font=dict(size=13, color='#ddd')),
        xaxis=dict(title='', color='#aaa', gridcolor='#333', automargin=True,
                   type='date', tickformat='%Y', tickformatstops=[
                       dict(dtickrange=[None, 'M3'],  value='%b %Y'),
                       dict(dtickrange=['M3', 'M18'], value='%b %Y'),
                       dict(dtickrange=['M18', None], value='%Y'),
                   ]),
        yaxis=dict(title='Depth (m)', autorange='reversed', color='#aaa', gridcolor='#333'),
        paper_bgcolor='#1a1a1a', plot_bgcolor='#111',
        margin=dict(l=60, r=20, t=50, b=60), autosize=True, dragmode='pan',
    )
    return fig


def make_source_table(bin_data):
    rows = []
    for ym, entry in sorted(bin_data['months'].items()):
        for s in entry['sources']:
            d = str(s['date'])
            rows.append({'Year-Month': f"{ym[0]}-{ym[1]:02d}", 'Date': f"{d[:4]}-{d[4:6]}-{d[6:8]}",
                         'Source': s.get('source', ''), 'Qual': s.get('qual', ''),
                         'Lon': f"{s['lon']:.3f}", 'Lat': f"{s['lat']:.3f}",
                         'File': s['file'], 'Prof index': s['prof_idx']})
    return rows


def _empty_fig(message=''):
    fig = go.Figure()
    fig.update_layout(
        paper_bgcolor='#1a1a1a', plot_bgcolor='#111',
        xaxis=dict(visible=False), yaxis=dict(visible=False),
        annotations=[dict(text=message, xref='paper', yref='paper', x=0.5, y=0.5,
                          showarrow=False, font=dict(size=14, color='#aaaaaa'))],
        autosize=True, margin=dict(l=20, r=20, t=20, b=20),
    )
    return fig

# ==============================================================================
# Layout
# ==============================================================================

app = Dash(__name__)
app.title = 'Profile Map'

DARK = '#1a1a1a'; MID = '#252525'; TEXT = '#dddddd'; ACC = '#88ccff'; DIM = '#aaaaaa'

YEAR_MIN = DATASET_FIRST_YM[0]
YEAR_MAX = DATASET_LAST_YM[0]
_YEAR_MARKS = {y: {'label': str(y), 'style': {'color': '#aaa', 'fontSize': '10px'}}
               for y in range(YEAR_MIN, YEAR_MAX + 1, 5)}
_SOURCE_OPTIONS = ([{'label': ' All (combined)', 'value': SOURCE_ALL}]
                   + [{'label': f' {s}', 'value': s} for s in SOURCES])
# every level is a tick (dots); depths written under every 16th, plus the last
_DEPTH_MARKS = {i: {'label': f'{depth[i]:.0f}', 'style': {'color': '#aaa', 'fontSize': '9px'}}
                for i in sorted(set(range(0, N_DEPTH, 16)) | {N_DEPTH - 1})}

_LABEL = {'fontSize': '12px', 'color': DIM, 'marginBottom': '6px', 'display': 'block'}
_STEP_BTN = {'fontSize': '11px', 'padding': '2px 5px', 'cursor': 'pointer', 'backgroundColor': '#333',
             'color': '#ddd', 'border': '1px solid #555', 'borderRadius': '3px', 'flex': '0 0 auto'}
_RADIO_BLOCK = dict(labelStyle={'display': 'block', 'fontSize': '12px', 'marginBottom': '5px',
                                'cursor': 'pointer', 'color': TEXT},
                    inputStyle={'marginRight': '6px'})
_RADIO_INLINE = dict(labelStyle={'display': 'inline-block', 'fontSize': '12px', 'marginRight': '12px',
                                 'cursor': 'pointer', 'color': TEXT},
                     inputStyle={'marginRight': '4px'})


def _control(label, child, cid=None, show=True):
    kw = {'id': cid} if cid else {}
    return html.Div([html.Label(label, style=_LABEL)] + (child if isinstance(child, list) else [child]),
                    style={'display': 'block' if show else 'none'}, **kw)


def _zt_graph(gid, message=''):
    return dcc.Loading(type='circle', color=ACC, children=[
        dcc.Graph(id=gid, figure=_empty_fig(message),
                  config={'scrollZoom': True, 'displayModeBar': True,
                          'modeBarButtonsToAdd': ['resetScale2d'],
                          'modeBarButtonsToRemove': ['toImage', 'sendDataToCloud'],
                          'displaylogo': False},
                  style={'height': '100%'})])


_INITIAL_PROJ = mc.PROJ_WORLD

app.layout = html.Div(
    style={'backgroundColor': DARK, 'color': TEXT, 'fontFamily': 'sans-serif', 'height': '100vh',
           'display': 'flex', 'flexDirection': 'column', 'padding': '8px', 'gap': '8px'},
    children=[

    # Header
    html.Div(style={'backgroundColor': MID, 'padding': '8px 14px', 'borderRadius': '4px',
                    'display': 'flex', 'alignItems': 'center', 'gap': '20px'}, children=[
        html.H2('Profile Map — Interactive Bin Explorer',
                style={'margin': 0, 'fontSize': '16px', 'color': ACC}),
        html.Span(('Pick a projection and view, then click a bin or a profile.  '
                   + (f'Data: {DATASET_FIRST_YM[0]}-{DATASET_FIRST_YM[1]:02d} – '
                      f'{DATASET_LAST_YM[0]}-{DATASET_LAST_YM[1]:02d}  |  '
                      if DATASET_FIRST_YM and DATASET_LAST_YM else '')
                   + f'{len(bins):,} bins  |  {", ".join(SOURCES)}'),
                  style={'fontSize': '12px', 'color': DIM}),
        html.Span(id='count-label', style={'fontSize': '12px', 'color': ACC, 'marginLeft': 'auto'}),
    ]),

    # Map row
    html.Div(style={'display': 'flex', 'gap': '8px', 'flex': '0 0 45vh', 'minHeight': '0'}, children=[

        html.Div(style={'flex': '1', 'backgroundColor': MID, 'borderRadius': '4px', 'padding': '4px',
                        'minHeight': '0', 'position': 'relative'}, children=[
            dcc.Graph(id='map', figure=make_map(_INITIAL_PROJ, SHOW_BINS, RES_FINE, SOURCE_ALL,
                                                f'{_INITIAL_PROJ}-0'),
                      config={'scrollZoom': True, 'doubleClick': False,
                              'displaylogo': False, 'displayModeBar': False},
                      style={'height': '100%'}),
            html.Button('⌂ Reset view', id='home-btn',
                        title='Back to the opening zoom and orientation (keeps depth, source and the rest)',
                        style={'position': 'absolute', 'top': '12px', 'right': '12px', 'zIndex': 10,
                               'fontSize': '13px', 'lineHeight': '1', 'padding': '6px 10px',
                               'cursor': 'pointer', 'backgroundColor': 'rgba(20,20,20,0.8)',
                               'color': '#ccc', 'border': '1px solid #666', 'borderRadius': '4px'}),
        ]),

        # Controls
        html.Div(style={'width': '240px', 'backgroundColor': MID, 'borderRadius': '4px',
                        'padding': '14px', 'display': 'flex', 'flexDirection': 'column',
                        'gap': '14px', 'overflowY': 'auto', 'minHeight': '0'}, children=[
            _control('Projection', dcc.RadioItems(
                id='proj-radio', value=_INITIAL_PROJ, **_RADIO_BLOCK,
                options=[{'label': f" {p['label']}", 'value': k} for k, p in mc.PROJECTIONS.items()])),
            _control('Show', dcc.RadioItems(
                id='show-radio', value=SHOW_BINS, **_RADIO_INLINE,
                options=[{'label': ' Bins', 'value': SHOW_BINS},
                         {'label': ' Profiles', 'value': SHOW_PROFILES}])),
            _control('Grid resolution', dcc.RadioItems(
                id='res-radio', value=RES_FINE, **_RADIO_BLOCK,
                options=[{'label': ' Fine (10242 bins)', 'value': RES_FINE},
                         {'label': ' Coarse (2562 bins)', 'value': RES_COARSE}])),
            _control('Source', dcc.Dropdown(id='source-select', options=_SOURCE_OPTIONS,
                                            value=SOURCE_ALL, clearable=False,
                                            style={'fontSize': '12px', 'color': '#111'})),

            # Bins view only
            html.Div(id='bins-controls', style={'display': 'flex', 'flexDirection': 'column',
                                                 'gap': '14px'}, children=[
                _control('Color by', dcc.Dropdown(id='cmode-select', options=CMODE_OPTIONS,
                                                  value=CMODE_PROB, clearable=False,
                                                  style={'fontSize': '12px', 'color': '#111'})),
                _control('Depth level', [
                    html.Div(id='depth-label', style={'fontSize': '11px', 'color': ACC,
                                                      'marginBottom': '4px'}),
                    html.Div(style={'display': 'flex', 'alignItems': 'center', 'gap': '2px'}, children=[
                        html.Button('◀', id='depth-up', title='One level shallower', style=_STEP_BTN),
                        html.Div(dcc.Slider(id='depth-slider', min=0, max=N_DEPTH - 1, step=1, value=0, dots=True,
                                            included=False, marks=_DEPTH_MARKS, updatemode='mouseup'),
                                 style={'flex': '1', 'minWidth': '0'}),
                        html.Button('▶', id='depth-down', title='One level deeper', style=_STEP_BTN),
                    ]),
                    html.Div('depth (m)', style={'fontSize': '9px', 'color': DIM, 'textAlign': 'right'}),
                    dcc.Checklist(id='rescale-depth', value=[],
                                  options=[{'label': ' Rescale colors to this depth', 'value': 'on'}],
                                  labelStyle={'fontSize': '12px', 'cursor': 'pointer', 'color': TEXT},
                                  inputStyle={'marginRight': '6px'}, style={'marginTop': '6px'}),
                    html.Div(f'(switches off at {RESCALE_RESET_DEPTH:.0f} m and above, and with Reset view)',
                             style={'fontSize': '10px', 'color': DIM}),
                ], cid='depth-slider-container', show=False),
                dcc.Checklist(id='edge-std', value=['on'],
                              options=[{'label': ' Edges: anomaly std', 'value': 'on'}],
                              labelStyle={'fontSize': '12px', 'cursor': 'pointer', 'color': TEXT},
                              inputStyle={'marginRight': '6px'}),
                _control('Probability: a month counts if the bin has', dcc.RadioItems(
                    id='prob-radio', value='>=1', **_RADIO_BLOCK,
                    options=[{'label': f' at least {v} profile' + ('s' if v > 1 else ''), 'value': k}
                             for k, v in PROB_KEYS.items()])),
                _control('Probability: out of', dcc.RadioItems(
                    id='basis-radio', value=BASIS_DAY0, **_RADIO_BLOCK,
                    options=[{'label': f' all {DATASET_MONTHS} months of the record', 'value': BASIS_DAY0},
                             {'label': " the months since the bin's first profile", 'value': BASIS_SINCE_FIRST}])),
                html.Div(id='prob-help', style={'fontSize': '11px', 'color': DIM, 'lineHeight': '1.45'}),
            ]),

            # Profiles view only
            _control('Time window', [
                html.Div(id='year-range-label', style={'fontSize': '11px', 'color': ACC,
                                                       'marginBottom': '4px'}),
                dcc.RangeSlider(id='year-slider', min=YEAR_MIN, max=YEAR_MAX, step=1,
                                value=[YEAR_MIN, YEAR_MAX], marks=_YEAR_MARKS, allowCross=False,
                                tooltip={'placement': 'bottom', 'always_visible': False}),
            ], cid='profiles-controls', show=False),

            # Bin panels
            _control('Field', dcc.RadioItems(
                id='field-radio', value='T', **_RADIO_BLOCK,
                options=[{'label': ' Temperature (T)', 'value': 'T'},
                         {'label': ' Salinity (S)', 'value': 'S'}])),
            _control('Anomaly reference', dcc.RadioItems(
                id='anom-radio', value='minus_clim', **_RADIO_BLOCK,
                options=[{'label': ' − climatology', 'value': 'minus_clim'},
                         {'label': ' − bin time-mean', 'value': 'anom'}])),
            dcc.Checklist(id='show-sources', value=[],
                          options=[{'label': ' Show source table', 'value': 'show'}],
                          labelStyle={'fontSize': '12px', 'cursor': 'pointer', 'color': TEXT},
                          inputStyle={'marginRight': '6px'}),
            html.Div(id='bin-info', style={'fontSize': '11px', 'color': DIM, 'lineHeight': '1.6'}),
        ]),
    ]),

    # Bin panels
    html.Div(style={'display': 'flex', 'gap': '8px', 'flex': '1', 'minHeight': '0',
                    'overflow': 'hidden'}, children=[
        html.Div(style={'flex': '1', 'backgroundColor': MID, 'borderRadius': '4px', 'padding': '4px',
                        'minWidth': '0'}, children=[_zt_graph('zt-obs', 'Click a bin or a profile on the map')]),
        html.Div(id='zt-anom-container', style={'flex': '1', 'backgroundColor': MID, 'borderRadius': '4px',
                                                'padding': '4px', 'minWidth': '0'},
                 children=[_zt_graph('zt-anom')]),
        html.Div(id='source-container',
                 style={'display': 'none', 'flex': '1', 'backgroundColor': MID, 'borderRadius': '4px',
                        'padding': '8px', 'flexDirection': 'column', 'overflow': 'hidden', 'minWidth': '0'},
                 children=[
            html.Div('Source profiles', style={'fontSize': '12px', 'color': ACC, 'marginBottom': '6px'}),
            dcc.Loading(type='circle', color=ACC, children=[
                html.Div(id='source-table-container', style={'flex': '1', 'overflowY': 'auto'})]),
        ]),
    ]),

    dcc.Store(id='selected-bin',    data=None),
    dcc.Store(id='selected-lonlat', data=None),     # Profiles view: the clicked dot
    dcc.Store(id='grid-res',        data=RES_FINE),
    dcc.Store(id='map-state',       data=mc.default_view(_INITIAL_PROJ)),
    dcc.Store(id='map-built',       data=0),        # bumped on each full redraw
    dcc.Store(id='subpoly-sig',     data=None),     # what the sub-polygon layer last showed
    dcc.Store(id='dots-sig',        data=None),     # where the dots were last thinned for
])

# ==============================================================================
# Callbacks
# ==============================================================================

@app.callback(
    Output('bins-controls', 'style'),
    Output('profiles-controls', 'style'),
    Output('depth-slider-container', 'style'),
    Output('depth-label', 'children'),
    Output('year-range-label', 'children'),
    Output('prob-help', 'children'),
    Input('show-radio', 'value'),
    Input('cmode-select', 'value'),
    Input('depth-slider', 'value'),
    Input('year-slider', 'value'),
    Input('prob-radio', 'value'),
    Input('basis-radio', 'value'),
    Input('source-select', 'value'),
)
def toggle_controls(show, cmode, depth_idx, years, prob_label, basis, source_sel):
    bins_view = show != SHOW_PROFILES
    depth_shown = bins_view and cmode != CMODE_PROB
    y0, y1 = years or [YEAR_MIN, YEAR_MAX]
    return ({'display': 'flex' if bins_view else 'none', 'flexDirection': 'column', 'gap': '14px'},
            {'display': 'none' if bins_view else 'block'},
            {'display': 'block' if depth_shown else 'none'},
            (f'level {(depth_idx or 0) + 1} of {N_DEPTH}: {depth[depth_idx or 0]:.0f} m'
             if depth_shown else ''),
            f'{y0} – {y1}',
            _prob_help(PROB_KEYS.get(prob_label, 1), basis, source_sel, cmode))


def _prob_help(n, basis, source_sel, cmode):
    """Plain-language definition of the probability, for the current choices."""
    noun = 'profile' if n == 1 else 'profiles'
    who = noun if source_sel == SOURCE_ALL else f'{source_sel} {noun}'
    out_of = (f'all {DATASET_MONTHS} months of the record '
              f'({DATASET_FIRST_YM[0]}-{DATASET_FIRST_YM[1]:02d} to {DATASET_LAST_YM[0]}-{DATASET_LAST_YM[1]:02d})'
              if basis != BASIS_SINCE_FIRST else 'the months from its first profile to the end of the record')
    text = (f'P(≥{n}) = the fraction of months in which a bin had at least {n} {who} '
            f'(at any depth), out of {out_of}. P = 0.8: 80% of those months had data.')
    return text if cmode == CMODE_PROB else text + ' Shown in the hover and bin info.'


@app.callback(
    Output('depth-slider', 'value'),
    Input('depth-up',   'n_clicks'),
    Input('depth-down', 'n_clicks'),
    State('depth-slider', 'value'),
    prevent_initial_call=True,
)
def step_depth(_up, _down, depth_idx):
    """◀ one level shallower, ▶ one level deeper (stops at the ends)."""
    step = -1 if ctx.triggered_id == 'depth-up' else 1
    new = int(np.clip((depth_idx or 0) + step, 0, N_DEPTH - 1))
    return new if new != depth_idx else no_update


@app.callback(
    Output('rescale-depth', 'value'),
    Input('depth-slider', 'value'),
    Input('home-btn', 'n_clicks'),
    State('rescale-depth', 'value'),
    prevent_initial_call=True,
)
def reset_rescale(depth_idx, _home, value):
    """'Rescale colors to this depth' switches itself off with Reset view, and when the
    depth comes back up to RESCALE_RESET_DEPTH or shallower."""
    if 'on' not in (value or []):
        return no_update
    if ctx.triggered_id == 'home-btn' or float(depth[depth_idx or 0]) <= RESCALE_RESET_DEPTH:
        return []
    return no_update


@app.callback(
    Output('map',          'figure'),
    Output('map-state',    'data'),
    Output('grid-res',     'data'),
    Output('selected-bin', 'data', allow_duplicate=True),
    Output('count-label',  'children'),
    Output('map-built',    'data'),
    Input('proj-radio',    'value'),
    Input('show-radio',    'value'),
    Input('res-radio',     'value'),
    Input('home-btn',      'n_clicks'),
    Input('source-select', 'value'),
    Input('year-slider',   'value'),
    State('prob-radio',    'value'),
    State('basis-radio',   'value'),
    State('cmode-select',  'value'),
    State('depth-slider',  'value'),
    State('selected-bin',  'data'),
    State('selected-lonlat', 'data'),
    State('map-state',     'data'),
    State('grid-res',      'data'),
    State('edge-std',      'value'),
    State('map-built',     'data'),
    State('rescale-depth', 'value'),
    prevent_initial_call=True,
)
def rebuild_map(proj, show, res, home_clicks, source_sel, years, prob_label, basis, cmode, depth_idx,
                selected_bin, sel_lonlat, map_state, prev_res, edge_std, built, rescale_v):
    """Full redraw: projection, view, grid, home button; and, in Profiles view, the
    source or year window (the dots change). In Bins view a source change is a
    cheap patch instead (on_bin_controls)."""
    trig = ctx.triggered_id
    if show != SHOW_PROFILES and trig in ('source-select', 'year-slider'):
        return no_update, no_update, no_update, no_update, no_update, no_update
    res = res or RES_FINE
    # zoom is kept across view / grid / data changes; reset by projection and home
    state = map_state if trig not in ('proj-radio', 'home-btn') else mc.default_view(proj)
    if not state or state.get('proj') != proj:
        state = mc.default_view(proj)
    state = no_update if state is map_state else state
    if res != prev_res:
        selected_bin, sel_lonlat = None, None
    common = dict(n_threshold=PROB_KEYS.get(prob_label, 1), basis=basis or BASIS_DAY0,
                  cmode=cmode or CMODE_PROB, depth_idx=depth_idx or 0)
    uirev = f'{proj}-{home_clicks or 0}'
    if show == SHOW_PROFILES:
        y0, y1 = years or [YEAR_MIN, YEAR_MAX]
        dots = query_dots(source_sel, res, y0, y1)
        count = (f'showing {dots["n_shown"]:,} of {dots["n_total"]:,} profiles'
                 if dots else 'no profiles in this window')
        fig = make_map(proj, show, res, source_sel, uirev, sel_lonlat=sel_lonlat, dots=dots, **common)
    else:
        count = ''
        fig = make_map(proj, show, res, source_sel, uirev, selected_bin=selected_bin,
                       show_edges='on' in (edge_std or []),
                       rescale=('on' in (rescale_v or [])) and trig != 'home-btn',
                       state=map_state if state is no_update else state, **common)
    # The user's zoom / pan survives this redraw because uirevision is unchanged
    # (it changes only with the projection and the home button).
    return fig, state, res, selected_bin, count, (built or 0) + 1


def _patch_edges(patched, ids, cmode, depth_idx, res, edge_on, bucket, uniform=False, src=SOURCE_ALL,
                 rescale=False):
    """Outline colors / widths (trace 4), their color bar (5) and the note (annotation 0)."""
    e = _edges(ids, cmode, depth_idx, res, bucket, uniform, src, rescale) if edge_on else None
    patched['data'][4]['visible'] = e is not None
    patched['data'][5]['visible'] = bool(e and e['cmax'] is not None)
    if e is not None:                       # geometry unchanged: only ids, colors, widths
        patched['data'][4]['locations'] = ids
        patched['data'][4]['z'] = [0] * len(ids)
        patched['data'][4]['marker']['line']['color'] = e['colors']
        patched['data'][4]['marker']['line']['width'] = e['widths']
        if e['cmax'] is not None:
            patched['data'][5]['marker']['cmax'] = e['cmax']
            patched['data'][5]['marker']['colorbar']['title']['text'] = e['title']
    patched['layout']['annotations'][0]['text'] = _edge_note(e)


@app.callback(
    Output('map', 'figure', allow_duplicate=True),
    Input('source-select', 'value'),
    Input('prob-radio',    'value'),
    Input('basis-radio',   'value'),
    Input('cmode-select',  'value'),
    Input('depth-slider',  'value'),
    Input('rescale-depth', 'value'),
    State('show-radio',    'value'),
    State('grid-res',      'data'),
    prevent_initial_call=True,
)
def on_bin_controls(source_sel, prob_label, basis, cmode, depth_idx, rescale_v, show, res):
    """Bins view: recolor the bins (trace 0) in place. (The outlines belong to
    on_subpolygons, which sees the zoom.)"""
    if show == SHOW_PROFILES:
        return no_update
    res = res or RES_FINE
    a = _bin_arrays(source_sel, PROB_KEYS[prob_label], basis, cmode or CMODE_PROB,
                    depth_idx or 0, res, 'on' in (rescale_v or []))
    patched = Patch()
    for key in ('z', 'colorscale', 'zmin', 'zmax'):
        patched['data'][0][key] = a[key]
    patched['data'][0]['locations']  = a['ids']
    patched['data'][0]['customdata'] = a['ids']
    patched['data'][0]['text']       = a['hover']
    patched['data'][0]['colorbar']['title']['text'] = a['cb_title']
    return patched


@app.callback(
    Output('map', 'figure', allow_duplicate=True),
    Output('count-label', 'children', allow_duplicate=True),
    Output('subpoly-sig', 'data'),
    Input('map-state',     'data'),
    Input('map-built',     'data'),
    Input('cmode-select',  'value'),
    Input('depth-slider',  'value'),
    Input('source-select', 'value'),
    Input('edge-std',      'value'),
    Input('rescale-depth', 'value'),
    State('show-radio',    'value'),
    State('grid-res',      'data'),
    State('subpoly-sig',   'data'),
    prevent_initial_call=True,
)
def on_subpolygons(state, built, cmode, depth_idx, source_sel, edge_std, rescale_v, show, res, prev_sig):
    """Bins view: fill / clear the sub-polygon trace (3) for the bins in view, and set
    the outlines (4, 5): density widths stepping with zoom, or uniform while sub-polygons
    show. Sends nothing when neither changed, so panning and zooming don't make Plotly
    redraw mid-gesture (that made the map bounce)."""
    if show == SHOW_PROFILES or not state:
        return no_update, no_update, no_update
    res, cmode, depth_idx = res or RES_FINE, cmode or CMODE_PROB, depth_idx or 0
    rescale = 'on' in (rescale_v or [])
    layer = _subpoly_layer(state, source_sel, cmode, depth_idx, res, rescale)
    gj, locs, z, hover, zmin, zmax, cs, note = layer
    bucket = mc.zoom_bucket(state)
    edge_on = 'on' in (edge_std or [])
    sig = [built, res, cmode, depth_idx, source_sel, bucket, len(locs), locs[0] if locs else '',
           locs[-1] if locs else '', note, edge_on, rescale]
    if sig == prev_sig:
        return no_update, no_update, no_update
    patched = Patch()
    patched['data'][3]['geojson'] = gj
    patched['data'][3]['locations'] = locs
    patched['data'][3]['z'] = z
    patched['data'][3]['text'] = hover
    patched['data'][3]['customdata'] = [int(f.split(':')[0]) for f in locs]
    patched['data'][3]['colorscale'] = cs
    patched['data'][3]['zmin'] = zmin
    patched['data'][3]['zmax'] = zmax
    # while the cells show, their color bar replaces the bins' (which they cover)
    patched['data'][3]['showscale'] = bool(locs)
    patched['data'][3]['colorbar']['title']['text'] = _subpoly_cbar_title(cmode)
    patched['data'][0]['showscale'] = not locs
    _patch_edges(patched, _visible_ids(source_sel, res), cmode, depth_idx, res, edge_on, bucket,
                 uniform=bool(locs), src=source_sel, rescale=rescale)
    return patched, note, sig


@app.callback(
    Output('map', 'figure', allow_duplicate=True),
    Output('count-label', 'children', allow_duplicate=True),
    Output('dots-sig', 'data'),
    Input('map-state',     'data'),
    Input('map-built',     'data'),
    State('show-radio',    'value'),
    State('source-select', 'value'),
    State('year-slider',   'value'),
    State('grid-res',      'data'),
    State('dots-sig',      'data'),
    prevent_initial_call=True,
)
def on_dots_view(state, built, show, source_sel, years, res, prev_sig):
    """Profiles view: re-thin the dots for the area in view, so zooming in shows up to
    TARGET_N of the profiles there (not just those kept from the whole-map thinning).
    Only redrawn when the view moves by more than a fraction of its size."""
    if show != SHOW_PROFILES or not state:
        return no_update, no_update, no_update
    z = mc.zoom_factor(state)
    step = 45.0 / z                                   # recenter grid: 1/4 of the view height
    sig = [built, state['proj'], round(np.log2(z) * 2) / 2,
           round(state['lon'] / step), round(state['lat'] / step)]
    if sig == prev_sig:
        return no_update, no_update, no_update
    if ctx.triggered_id == 'map-built' and z <= 1.01:   # the redraw already thinned for the whole map
        return no_update, no_update, sig
    y0, y1 = years or [YEAR_MIN, YEAR_MAX]
    dots = query_dots(source_sel, res or RES_FINE, y0, y1, view=state if z > 1.01 else None)
    patched = Patch()
    for i, tr in enumerate(_dot_traces(dots, source_sel)):
        for key in ('lat', 'lon', 'text', 'customdata'):
            patched['data'][1 + i][key] = tr[key] if tr[key] is not None else []
    count = (f'showing {dots["n_shown"]:,} of {dots["n_total"]:,} profiles'
             + (' in view' if z > 1.01 else '') if dots else 'no profiles in view')
    return patched, count, sig


@app.callback(
    Output('selected-bin',    'data'),
    Output('selected-lonlat', 'data'),
    Output('map', 'figure', allow_duplicate=True),
    Input('map', 'clickData'),
    State('selected-bin',  'data'),
    State('show-radio',    'value'),
    State('source-select', 'value'),
    State('grid-res',      'data'),
    prevent_initial_call=True,
)
def on_map_click(click_data, current_bin, show, source_sel, res):
    points = (click_data or {}).get('points', [])
    if not points:
        return no_update, no_update, no_update
    p = points[0]
    bid = p.get('customdata', p.get('location'))
    res = res or RES_FINE
    if bid is None or bid not in _grid(res)[0]:
        return no_update, no_update, no_update
    patched = Patch()
    if show == SHOW_PROFILES:
        sel = {'lon': p.get('lon'), 'lat': p.get('lat')}
        i = 1 + _n_dot_traces(source_sel or SOURCE_ALL)       # the star (see profiles_traces)
        patched['data'][i]['lon'] = [sel['lon']]
        patched['data'][i]['lat'] = [sel['lat']]
        patched['data'][i]['marker'] = dict(size=14, symbol='star', color='white',
                                            line=dict(color='black', width=1))
        return bid, sel, patched
    lon, lat = _ring_lonlat(bid, res)
    patched['data'][1]['lon'] = lon
    patched['data'][1]['lat'] = lat
    return bid, no_update, patched


@app.callback(
    Output('map', 'figure', allow_duplicate=True),
    Output('map-state', 'data', allow_duplicate=True),
    Input('map', 'relayoutData'),
    State('map-state', 'data'),
    State('proj-radio', 'value'),
    prevent_initial_call=True,
)
def on_map_view(relayout, state, proj):
    """Track zoom / pan, and push it back inside the allowed range."""
    new_state, clamp = mc.update_view(relayout, state, proj)
    if not clamp:
        return no_update, new_state
    patched = Patch()
    if 'scale' in clamp:
        patched['layout']['geo']['projection']['scale'] = clamp['scale']
    if 'center_lat' in clamp:
        patched['layout']['geo']['center']['lat'] = clamp['center_lat']
    if 'rotation_lat' in clamp:
        patched['layout']['geo']['projection']['rotation']['lat'] = clamp['rotation_lat']
    return patched, new_state


@app.callback(
    Output('zt-anom-container', 'style'),
    Output('source-container',  'style'),
    Input('show-sources', 'value'),
)
def toggle_right_panel(show_sources):
    show = 'show' in (show_sources or [])
    anom_style = {'flex': '1', 'backgroundColor': MID, 'borderRadius': '4px', 'padding': '4px',
                  'minWidth': '0', 'display': 'none' if show else 'block'}
    src_style  = {'flex': '1', 'backgroundColor': MID, 'borderRadius': '4px', 'padding': '8px',
                  'minWidth': '0', 'overflow': 'hidden', 'flexDirection': 'column',
                  'display': 'flex' if show else 'none'}
    return anom_style, src_style


def _axis_patch_from_relayout(relayout):
    if not relayout:
        return None
    patch = Patch(); touched = False
    if 'xaxis.range[0]' in relayout and 'xaxis.range[1]' in relayout:
        patch['layout']['xaxis']['range'] = [relayout['xaxis.range[0]'], relayout['xaxis.range[1]']]
        patch['layout']['xaxis']['autorange'] = False; touched = True
    if 'yaxis.range[0]' in relayout and 'yaxis.range[1]' in relayout:
        patch['layout']['yaxis']['range'] = [relayout['yaxis.range[0]'], relayout['yaxis.range[1]']]
        patch['layout']['yaxis']['autorange'] = False; touched = True
    if relayout.get('xaxis.autorange'):
        patch['layout']['xaxis']['autorange'] = True; touched = True
    if relayout.get('yaxis.autorange'):
        patch['layout']['yaxis']['autorange'] = 'reversed'; touched = True
    return patch if touched else None


@app.callback(Output('zt-anom', 'figure', allow_duplicate=True), Input('zt-obs', 'relayoutData'),
              prevent_initial_call=True)
def sync_obs_to_anom(relayout):
    p = _axis_patch_from_relayout(relayout)
    return p if p is not None else no_update


@app.callback(Output('zt-obs', 'figure', allow_duplicate=True), Input('zt-anom', 'relayoutData'),
              prevent_initial_call=True)
def sync_anom_to_obs(relayout):
    p = _axis_patch_from_relayout(relayout)
    return p if p is not None else no_update


@app.callback(
    Output('zt-obs',                 'figure'),
    Output('zt-anom',                'figure'),
    Output('source-table-container', 'children'),
    Output('bin-info',               'children'),
    Input('selected-bin',  'data'),
    Input('source-select', 'value'),
    Input('field-radio',   'value'),
    Input('anom-radio',    'value'),
    Input('prob-radio',    'value'),
    Input('basis-radio',   'value'),
    Input('grid-res',      'data'),
)
def update_plots(bin_id, source_sel, field, anom_ref, prob_label, basis, res):
    res = res or RES_FINE
    b_dict, c_dict, _, _ = _grid(res)
    res_label = '10242-bin' if res == RES_FINE else '2562-bin'

    if bin_id is None or bin_id not in b_dict:
        return (_empty_fig('Click a bin or a profile on the map'), _empty_fig(''),
                html.Div('No bin selected.', style={'color': DIM, 'fontSize': '12px'}), '')

    resolved = resolve_bin(b_dict[bin_id], source_sel, res)
    lon, lat = c_dict.get(bin_id, (np.nan, np.nan))
    n_months = len(resolved['months'])
    src_note = 'All sources' if source_sel == SOURCE_ALL else source_sel
    if n_months == 0:
        msg = _empty_fig(f'Bin {bin_id} ({res_label}) has no {src_note.lower()} data')
        return (msg, _empty_fig(''),
                html.Div(f'No {src_note} profiles in this bin.', style={'color': DIM, 'fontSize': '12px'}),
                html.Div(f'Bin {bin_id}: no {src_note} data', style={'color': DIM}))

    obs_fig  = make_zt_figure(resolved, field, depth)
    anom_fig = make_zt_figure(resolved, f'{field}_minus_clim' if anom_ref == 'minus_clim' else f'{field}_anom',
                              depth)
    table = dash_table.DataTable(
        data=make_source_table(resolved),
        columns=[{'name': c, 'id': c} for c in
                 ('Year-Month', 'Date', 'Source', 'Qual', 'Lon', 'Lat', 'File', 'Prof index')],
        style_table={'overflowY': 'auto', 'height': '100%'},
        style_cell={'backgroundColor': '#1e1e1e', 'color': '#ccc', 'fontSize': '11px',
                    'padding': '3px 7px', 'border': '1px solid #333', 'textAlign': 'left',
                    'whiteSpace': 'nowrap', 'overflow': 'hidden', 'textOverflow': 'ellipsis',
                    'maxWidth': '180px'},
        style_header={'backgroundColor': '#2d2d2d', 'color': ACC, 'fontWeight': 'bold', 'fontSize': '11px'},
        style_data_conditional=[{'if': {'row_index': 'odd'}, 'backgroundColor': '#232323'}],
        page_size=200, sort_action='native', filter_action='native',
    )
    n_threshold = PROB_KEYS[prob_label]
    prob_val    = bin_probability(b_dict[bin_id], source_sel, n_threshold, basis)
    info = [
        html.Div(f'Bin {bin_id} ({res_label})', style={'color': ACC, 'fontWeight': 'bold'}),
        html.Div(src_note, style={'color': DIM}),
        html.Div(f'lon {lon:.2f}°'), html.Div(f'lat {lat:.2f}°'),
        html.Div(f'{n_months} months'), html.Div(f'{resolved["n_profiles"]} profiles'),
        html.Div(f'{len(b_dict[bin_id]["by_source"])} source(s) in bin', style={'color': DIM}),
        html.Div(f'P(≥{n_threshold}) = {prob_val:.3f}'),
        html.Div('(since first obs)' if basis == BASIS_SINCE_FIRST else '(full span)', style={'color': DIM}),
    ]
    return obs_fig, anom_fig, table, info


if __name__ == '__main__':
    print("Starting app — open http://127.0.0.1:8050")
    app.run(debug=False, host='127.0.0.1', port=8050)
