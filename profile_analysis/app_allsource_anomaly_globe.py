"""
Multi-source interactive globe app with anomaly coloring.

Extends app_allsource_globe.py with two additional globe/flat-map color modes:
  - Anomaly mean  : long-run mean of (T or S) - climatology across all profiles in bin
  - Anomaly std   : std of the same, indicating temporal variability

A "Color by" selector switches between:
  - Visit probability  (original mode, no depth slider)
  - T anomaly mean     (depth slider appears)
  - T anomaly std
  - S anomaly mean
  - S anomaly std

Requires the data store built by build_allsource_bin_timeseries_with_anomalies.py,
which stores T_anom_mean / T_anom_std / S_anom_mean / S_anom_std per bin in zarr.

Run with:
    python app_allsource_anomaly_globe.py
Then open http://127.0.0.1:8050
"""

import json
import numpy as np
import plotly.graph_objects as go
import cartopy.feature as cfeature
import argparse
import sqlite3
import zarr
from dash import Dash, dcc, html, Input, Output, State, dash_table, Patch, no_update
from pathlib import Path
from scipy.spatial import SphericalVoronoi

# ==============================================================================
# Config
# ==============================================================================

_OUTPUT_BASE = Path('/Users/brucel/ecco/yip/profile_data/z_profile_file_analysis/'
                    'allsource_bin_timeseries_with_anomalies')

def _latest_store(base):
    """Return the most recent timestamped store directory, falling back to base itself."""
    parent = base.parent
    candidates = sorted(
        p for p in parent.glob(base.name + '_????????_??????')
        if p.is_dir()
    )
    return candidates[-1] if candidates else base

_app_parser = argparse.ArgumentParser(description='Allsource anomaly globe app.')
_app_parser.add_argument('--data', default=None,
                         help='Path to the store directory to load. '
                              'Defaults to the most recently timestamped '
                              'allsource_bin_timeseries_with_anomalies_* directory.')
_app_args, _unknown = _app_parser.parse_known_args()

_DATA_DIR     = Path(_app_args.data) if _app_args.data else _latest_store(_OUTPUT_BASE)
ZARR_STORE    = _DATA_DIR / (_DATA_DIR.name + '.zarr')
METADATA_JSON = _DATA_DIR / (_DATA_DIR.name + '_metadata.json')
print(f"Using data store: {_DATA_DIR}", flush=True)
GEODESIC_2562_CSV = Path('/Users/brucel/ecco/yip/sample_data/ecco-insitu/sweet_gdrive/geodesic/02562_bin_locations.csv')

SOURCE_ALL = 'ALL'

CELL_WINDING = -1

# Color-by modes
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

# Colorscale + label per mode
_CMODE_CS = {
    CMODE_PROB:   ('RdYlGn',  'P(≥1)',   0.0,  1.0,  False),  # (cscale, label, zmin, zmax, diverging)
    CMODE_T_MEAN: ('RdBu_r',  'T anom mean (°C)',  None, None, True),
    CMODE_T_STD:  ('cividis', 'T anom std (°C)',   0.0,  None, False),
    CMODE_S_MEAN: ('RdBu_r',  'S anom mean (PSU)', None, None, True),
    CMODE_S_STD:  ('cividis', 'S anom std (PSU)',  0.0,  None, False),
}

# ==============================================================================
# Load metadata
# ==============================================================================

print("Loading metadata JSON...", flush=True)
with open(METADATA_JSON, 'r') as _f:
    _meta = json.load(_f)

depth               = np.array(_meta['depth_list'], dtype=np.float32)
SOURCES             = _meta['sources']
DATASET_FIRST_YM    = tuple(_meta['dataset_first_ym']) if _meta.get('dataset_first_ym') else None
DATASET_LAST_YM     = tuple(_meta['dataset_last_ym'])  if _meta.get('dataset_last_ym')  else None
DATASET_MONTHS      = _meta.get('total_months', 1)
bin_centres         = {int(k): tuple(v) for k, v in _meta['bin_centres'].items()}
bin_centres_coarse  = {int(k): tuple(v) for k, v in _meta['bin_centres_coarse'].items()}
N_DEPTH             = len(depth)


def _parse_bin_meta(raw):
    months = {}
    for ym_str, m in raw['months'].items():
        y, mo = int(ym_str[:4]), int(ym_str[5:])
        by_src = {}
        for src, e in m['by_source'].items():
                by_src[src] = {'n': e['n']}
        months[(y, mo)] = {'n': m['n'], 'by_source': by_src}
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
print(f"  {len(bins):,} fine bins, {len(bins_coarse):,} coarse bins, "
      f"{len(SOURCES)} sources: {', '.join(SOURCES)}")

print("Opening zarr store...", flush=True)
_ZARR = zarr.open_group(str(ZARR_STORE), mode='r')
print("  Ready.")

# Open SQLite profile-records db (read-only; profile rows fetched lazily on bin click)
_PROFILES_DB = _DATA_DIR / (_DATA_DIR.name + '_profiles.db')
_DB_CON = sqlite3.connect(f'file:{_PROFILES_DB}?mode=ro', uri=True, check_same_thread=False)
_DB_CON.row_factory = sqlite3.Row


def _query_profiles(resolution, bid):
    """Return profile record dicts for a bin, grouped by {(year, month): {src: [recs]}}."""
    cur = _DB_CON.execute(
        'SELECT src, year, month, date, lon, lat, file, prof_idx, qual '
        'FROM profiles WHERE resolution=? AND bin_id=?',
        (resolution, bid)
    )
    result = {}
    for row in cur:
        ym = (row['year'], row['month'])
        src = row['src']
        result.setdefault(ym, {}).setdefault(src, []).append({
            'date':     row['date'],
            'lon':      row['lon'],
            'lat':      row['lat'],
            'file':     row['file'],
            'prof_idx': row['prof_idx'],
            'qual':     row['qual'],
            'source':   src,
        })
    return result

# ==============================================================================
# Zarr accessors
# ==============================================================================

def _zarr_src_monthly(resolution, bid, src):
    """Return (ym_index, arrays_dict) for one source.

    ym_index : (n_months,) int32  YYYYMM integers
    arrays_dict : key -> (n_months, N_DEPTH) float32
    """
    node     = _ZARR[f"{resolution}/{bid}/by_source/{src}"]
    ym_index = np.asarray(_ZARR[f"{resolution}/{bid}/ym_index"], dtype=np.int32)
    arrays   = {k: np.asarray(node[k], dtype=np.float32) for k in ('T', 'S', 'Tclim', 'Sclim')}
    return ym_index, arrays


def _zarr_arrays(resolution, bid, ym, src):
    """Read one month's per-source T/S/Tclim/Sclim. Returns dict of (N_DEPTH,) float32."""
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
    """Return (T_anom_mean, T_anom_std, S_anom_mean, S_anom_std) as float32 arrays."""
    node = _ZARR[f"{resolution}/{bid}"]
    return (np.asarray(node['T_anom_mean'], dtype=np.float32),
            np.asarray(node['T_anom_std'],  dtype=np.float32),
            np.asarray(node['S_anom_mean'], dtype=np.float32),
            np.asarray(node['S_anom_std'],  dtype=np.float32))

# ==============================================================================
# Geodesic cell polygons
# ==============================================================================

def _build_cell_geojson(centres):
    ids = sorted(centres.keys())
    lon = np.radians([centres[b][0] for b in ids])
    lat = np.radians([centres[b][1] for b in ids])
    xyz = np.stack([np.cos(lat) * np.cos(lon),
                    np.cos(lat) * np.sin(lon),
                    np.sin(lat)], axis=-1)
    sv = SphericalVoronoi(xyz, radius=1.0, center=np.array([0.0, 0.0, 0.0]))
    sv.sort_vertices_of_regions()

    def xyz2lonlat(v):
        return (np.degrees(np.arctan2(v[1], v[0])),
                np.degrees(np.arcsin(np.clip(v[2], -1.0, 1.0))))

    features = []
    for k, bid in enumerate(ids):
        region  = sv.regions[k]
        verts3d = sv.vertices[region]
        area_vec = np.zeros(3)
        for i in range(len(verts3d)):
            area_vec += np.cross(verts3d[i], verts3d[(i + 1) % len(verts3d)])
        if np.dot(area_vec, xyz[k]) * CELL_WINDING < 0:
            region = region[::-1]
        ring = [xyz2lonlat(sv.vertices[i]) for i in region]
        lons = np.array([p[0] for p in ring])
        lats = np.array([p[1] for p in ring])
        c_lon = centres[bid][0]
        lons  = c_lon + ((lons - c_lon + 180.0) % 360.0 - 180.0)
        coords = [[float(a), float(b)] for a, b in zip(lons, lats)]
        coords.append(coords[0])
        features.append({'type': 'Feature', 'id': int(bid),
                         'geometry': {'type': 'Polygon', 'coordinates': [coords]},
                         'properties': {}})
    return {'type': 'FeatureCollection', 'features': features}, ids


print("Building fine geodesic cell polygons...", flush=True)
CELL_GEOJSON, CELL_IDS = _build_cell_geojson(bin_centres)
print(f"  {len(CELL_IDS):,} fine cells built")

print("Building coarse geodesic cell polygons...", flush=True)
_coarse_raw = np.genfromtxt(GEODESIC_2562_CSV, delimiter=',')
_all_coarse_centres = {i + 1: (float(_coarse_raw[i, 0]), float(_coarse_raw[i, 1]))
                       for i in range(len(_coarse_raw))}
for k, v in _all_coarse_centres.items():
    bin_centres_coarse.setdefault(k, v)
COARSE_GEOJSON, COARSE_IDS = _build_cell_geojson(_all_coarse_centres)
print(f"  {len(COARSE_IDS):,} coarse cells built")

# ==============================================================================
# Helpers
# ==============================================================================

def _months_between(ym0, ym1):
    return (ym1[0] - ym0[0]) * 12 + (ym1[1] - ym0[1]) + 1

PROB_KEYS = {'>=1': 1, '>=2': 2, '>=3': 3}
BASIS_DAY0        = 'day0'
BASIS_SINCE_FIRST = 'since_first'
RES_FINE   = 'fine'
RES_COARSE = 'coarse'


def _grid(res):
    if res == RES_COARSE:
        return bins_coarse, bin_centres_coarse, COARSE_GEOJSON, COARSE_IDS
    return bins, bin_centres, CELL_GEOJSON, CELL_IDS


def _pool_level_weighted(entries, key):
    num = np.zeros(N_DEPTH)
    den = np.zeros(N_DEPTH)
    for e in entries:
        v = np.asarray(e[key], dtype=float)
        n = e['n']
        mask = ~np.isnan(v)
        num[mask] += v[mask] * n
        den[mask] += n
    out = np.full(N_DEPTH, np.nan)
    nz = den > 0
    out[nz] = num[nz] / den[nz]
    return out


def bin_has_source(bin_data, source_sel):
    return source_sel == SOURCE_ALL or source_sel in bin_data['by_source']


def bin_probability(bin_data, source_sel, n_threshold, basis):
    if source_sel == SOURCE_ALL:
        hits     = bin_data['hit_counts'].get(n_threshold, 0)
        first_ym = bin_data['first_ym']
    else:
        roll = bin_data['by_source'].get(source_sel)
        if roll is None:
            return 0.0
        hits     = roll['hit_counts'].get(n_threshold, 0)
        first_ym = roll['first_ym']
    if basis == BASIS_SINCE_FIRST and first_ym and DATASET_LAST_YM:
        denom = _months_between(first_ym, DATASET_LAST_YM)
    else:
        denom = DATASET_MONTHS
    return hits / denom if denom > 0 else 0.0


def _visible_ids(source_sel, res=RES_FINE):
    b_dict, c_dict, _, _ = _grid(res)
    return [b for b in sorted(b_dict.keys())
            if b in c_dict and bin_has_source(b_dict[b], source_sel)]


def resolve_bin(bin_data, source_sel, res=RES_FINE):
    bid  = bin_data['bin_id']
    zres = res

    # Fetch profile records for this bin from SQLite (one query for all months/sources)
    db_recs = _query_profiles(zres, bid)

    months = {}
    for ym, m in bin_data['months'].items():
        if source_sel != SOURCE_ALL and source_sel not in m['by_source']:
            continue
        src_list = (list(m['by_source'].keys()) if source_sel == SOURCE_ALL else [source_sel])
        entries  = []
        for src in src_list:
            if src not in m['by_source']:
                continue
            arrays = _zarr_arrays(zres, bid, ym, src)
            entry  = dict(arrays)
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
        T_mean, S_mean = _zarr_means(zres, bid)
        hit_counts  = bin_data['hit_counts']
        first_ym    = bin_data['first_ym']
        last_ym     = bin_data['last_ym']
        n_profiles  = sum(r['n_profiles'] for r in bin_data['by_source'].values())
    else:
        T_mean, S_mean = _zarr_means(zres, bid, source_sel)
        roll       = bin_data['by_source'][source_sel]
        hit_counts = roll['hit_counts']
        first_ym   = roll['first_ym']
        last_ym    = roll['last_ym']
        n_profiles = roll['n_profiles']

    return {
        'bin_id':     bin_data['bin_id'],
        'months':     months,
        'T_mean':     np.asarray(T_mean, dtype=float),
        'S_mean':     np.asarray(S_mean, dtype=float),
        'hit_counts': hit_counts,
        'first_ym':   first_ym,
        'last_ym':    last_ym,
        'n_profiles': n_profiles,
    }

# ==============================================================================
# Globe/flat map color arrays
# ==============================================================================

def _anom_arrays_for_display(ids_clean, cmode, depth_idx, res):
    """Return (z_values, colorscale, zmin, zmax, label, hover_suffix) for anomaly modes."""
    zres = res
    cs_info = _CMODE_CS[cmode]
    colorscale, label, zmin_fixed, zmax_fixed, diverging = cs_info

    z_vals = []
    for bid in ids_clean:
        try:
            T_am, T_as, S_am, S_as = _zarr_anom_stats(zres, bid)
            if cmode == CMODE_T_MEAN:
                v = float(T_am[depth_idx]) if not np.isnan(T_am[depth_idx]) else 0.0
            elif cmode == CMODE_T_STD:
                v = float(T_as[depth_idx]) if not np.isnan(T_as[depth_idx]) else 0.0
            elif cmode == CMODE_S_MEAN:
                v = float(S_am[depth_idx]) if not np.isnan(S_am[depth_idx]) else 0.0
            else:  # S_STD
                v = float(S_as[depth_idx]) if not np.isnan(S_as[depth_idx]) else 0.0
        except Exception:
            v = 0.0
        z_vals.append(v)

    z_arr = np.array(z_vals, dtype=float)

    if zmin_fixed is not None:
        zmin = zmin_fixed
    elif diverging:
        abs_max = np.nanpercentile(np.abs(z_arr[np.isfinite(z_arr)]), 97) if z_arr.size else 1.0
        zmin = -abs_max
    else:
        zmin = float(np.nanpercentile(z_arr[np.isfinite(z_arr)], 2)) if z_arr.size else 0.0

    if zmax_fixed is not None:
        zmax = zmax_fixed
    elif diverging:
        zmax = -zmin
    else:
        zmax = float(np.nanpercentile(z_arr[np.isfinite(z_arr)], 98)) if z_arr.size else 1.0

    depth_m = float(depth[depth_idx])
    hover_suffix = f'<br>{label} @ {depth_m:.0f} m = {{:.3f}}'
    return z_vals, colorscale, zmin, zmax, label, depth_m


def _globe_arrays(source_sel, n_threshold, basis, cmode, depth_idx, res=RES_FINE):
    b_dict, c_dict, _, _ = _grid(res)
    ids_clean  = _visible_ids(source_sel, res)
    prob_label = next(k for k, v in PROB_KEYS.items() if v == n_threshold)
    src_note   = 'All sources' if source_sel == SOURCE_ALL else source_sel

    if cmode == CMODE_PROB:
        probs = [bin_probability(b_dict[b], source_sel, n_threshold, basis) for b in ids_clean]
        hover = [
            f'Bin {b} — {src_note}'
            f'<br>lon={c_dict[b][0]:.2f}, lat={c_dict[b][1]:.2f}'
            f'<br>P({prob_label} prof/month) = {bin_probability(b_dict[b], source_sel, n_threshold, basis):.3f}'
            for b in ids_clean
        ]
        cs, zmin, zmax, cb_title = 'RdYlGn', 0.0, 1.0, f'P({prob_label})'
        return dict(ids=ids_clean, z=probs, hover=hover,
                    colorscale=cs, zmin=zmin, zmax=zmax, cb_title=cb_title)
    else:
        z_vals, colorscale, zmin, zmax, label, depth_m = _anom_arrays_for_display(
            ids_clean, cmode, depth_idx, res)
        probs = [bin_probability(b_dict[b], source_sel, n_threshold, basis) for b in ids_clean]
        hover = [
            f'Bin {b} — {src_note}'
            f'<br>lon={c_dict[b][0]:.2f}, lat={c_dict[b][1]:.2f}'
            f'<br>{label} @ {depth_m:.0f} m = {z:.3f}'
            f'<br>P({prob_label}) = {p:.3f}'
            for b, z, p in zip(ids_clean, z_vals, probs)
        ]
        return dict(ids=ids_clean, z=z_vals, hover=hover,
                    colorscale=colorscale, zmin=zmin, zmax=zmax, cb_title=f'{label}')


def _selected_trace(selected_bin_id, res=RES_FINE):
    if selected_bin_id is None:
        return go.Scattergeo(lon=[], lat=[], mode='lines', hoverinfo='skip',
                             showlegend=False, name='selection')
    _, _, geojson, _ = _grid(res)
    ring_lon, ring_lat = [], []
    for feat in geojson['features']:
        if feat['id'] == selected_bin_id:
            coords = feat['geometry']['coordinates'][0]
            ring_lon = [c[0] for c in coords]
            ring_lat = [c[1] for c in coords]
            break
    return go.Scattergeo(
        lon=ring_lon, lat=ring_lat, mode='lines',
        line=dict(color='rgba(255,255,255,0.95)', width=2.5),
        hoverinfo='skip', showlegend=False, name='selection',
    )


def _make_choropleth(a, geojson):
    return go.Choropleth(
        geojson=geojson,
        locations=a['ids'],
        z=a['z'],
        featureidkey='id',
        colorscale=a['colorscale'],
        zmin=a['zmin'], zmax=a['zmax'],
        marker=dict(line=dict(color='rgba(40,40,40,0.5)', width=0.3)),
        colorbar=dict(title=a['cb_title'], thickness=14, len=0.6),
        text=a['hover'],
        hoverinfo='text',
        customdata=a['ids'],
    )


GEO_STYLE = dict(
    showland=True,    landcolor='#2a2a2a',
    showocean=True,   oceancolor='#0d1b2a',
    showlakes=True,   lakecolor='#3a3a3a',
    showrivers=True,  rivercolor='#3a3a3a',
    showcoastlines=True, coastlinecolor='#cccccc', coastlinewidth=0.8,
    showframe=False,  bgcolor='#1a1a1a',
)


def make_globe(source_sel=SOURCE_ALL, selected_bin_id=None,
               n_threshold=1, basis=BASIS_DAY0, scale=1.0, res=RES_FINE,
               cmode=CMODE_PROB, depth_idx=0):
    _, _, geojson, _ = _grid(res)
    a   = _globe_arrays(source_sel, n_threshold, basis, cmode, depth_idx, res)
    fig = go.Figure(_make_choropleth(a, geojson))
    fig.add_trace(_selected_trace(selected_bin_id, res))
    fig.update_geos(projection_type='orthographic',
                    projection_scale=scale,
                    projection_rotation=dict(lon=0, lat=20, roll=0),
                    **GEO_STYLE)
    fig.update_layout(margin=dict(l=0, r=0, t=0, b=0),
                      paper_bgcolor='#1a1a1a', geo=dict(bgcolor='#1a1a1a'),
                      autosize=True, dragmode=False, uirevision='globe')
    return fig


def _build_coastline_trace():
    lons, lats = [], []
    for g in cfeature.NaturalEarthFeature('physical', 'coastline', '50m').geometries():
        parts = [g] if g.geom_type == 'LineString' else list(g.geoms)
        for part in parts:
            x, y = part.xy
            lons.extend(list(x) + [None]); lats.extend(list(y) + [None])
    return go.Scattergeo(lon=lons, lat=lats, mode='lines',
                         line=dict(color='#4a9eff', width=0.9),
                         hoverinfo='skip', showlegend=False)

_COAST_TRACE = _build_coastline_trace()


def make_flat(source_sel=SOURCE_ALL, selected_bin_id=None,
              n_threshold=1, basis=BASIS_DAY0, res=RES_FINE, uirevision='flat',
              cmode=CMODE_PROB, depth_idx=0):
    _, _, geojson, _ = _grid(res)
    a   = _globe_arrays(source_sel, n_threshold, basis, cmode, depth_idx, res)
    fig = go.Figure(_make_choropleth(a, geojson))
    fig.add_trace(_selected_trace(selected_bin_id, res))
    fig.add_trace(_COAST_TRACE)
    fig.update_geos(projection_type='natural earth',
                    showcoastlines=False,
                    lataxis_range=[-90, 90], lonaxis_range=[-180, 180],
                    **GEO_STYLE)
    fig.update_layout(margin=dict(l=0, r=0, t=0, b=0),
                      paper_bgcolor='#1a1a1a', geo=dict(bgcolor='#1a1a1a'),
                      autosize=True, dragmode='pan', uirevision=uirevision)
    return fig

# ==============================================================================
# Z-T heatmap
# ==============================================================================

def make_zt_figure(bin_data, var, depth):
    months_dict = bin_data['months']
    sorted_yms  = sorted(months_dict.keys())
    if not sorted_yms:
        return _empty_fig('No data for this source in this bin')
    x = [f'{ym[0]}-{ym[1]:02d}-15' for ym in sorted_yms]
    is_diff_clim = var in ('T_minus_clim', 'S_minus_clim')
    is_anom      = var in ('T_anom', 'S_anom')
    base_letter  = var[0]
    if is_diff_clim:
        obs_key, clim_key = base_letter, base_letter + 'clim'
        matrix = np.stack([months_dict[ym][obs_key] - months_dict[ym][clim_key]
                           for ym in sorted_yms], axis=0).T
        colorscale = 'RdBu_r'
        abs_max = np.nanpercentile(np.abs(matrix[np.isfinite(matrix)]), 97) if np.any(np.isfinite(matrix)) else 1
        zmin, zmax = -abs_max, abs_max
        title_var, units = f'{base_letter} − climatology', '°C' if base_letter == 'T' else 'PSU'
    elif is_anom:
        matrix = np.stack([months_dict[ym][base_letter] for ym in sorted_yms], axis=0).T
        matrix = matrix - bin_data[f'{base_letter}_mean'][:, np.newaxis]
        colorscale = 'RdBu_r'
        abs_max = np.nanpercentile(np.abs(matrix[np.isfinite(matrix)]), 97) if np.any(np.isfinite(matrix)) else 1
        zmin, zmax = -abs_max, abs_max
        title_var, units = f'{base_letter} − bin time-mean', '°C' if base_letter == 'T' else 'PSU'
    else:
        matrix = np.stack([months_dict[ym][var] for ym in sorted_yms], axis=0).T
        valid  = matrix[np.isfinite(matrix)]
        if base_letter == 'T':
            colorscale = 'Thermal'
            zmin = np.nanpercentile(valid, 2)  if len(valid) else -2
            zmax = np.nanpercentile(valid, 98) if len(valid) else 30
        else:
            colorscale = 'Haline'
            zmin = np.nanpercentile(valid, 2)  if len(valid) else 30
            zmax = np.nanpercentile(valid, 98) if len(valid) else 38
        clim_suffix = ' (clim)' if var.endswith('clim') else ''
        title_var, units = base_letter + clim_suffix, '°C' if base_letter == 'T' else 'PSU'

    fig = go.Figure(go.Heatmap(
        x=x, y=depth, z=matrix,
        colorscale=colorscale, zmin=zmin, zmax=zmax,
        colorbar=dict(title=units, thickness=12),
        hovertemplate='%{x|%b %Y}<br>Depth: %{y} m<br>Value: %{z:.3f}<extra></extra>',
    ))
    bid = bin_data['bin_id']
    lon, lat = (bin_centres.get(bid) or bin_centres_coarse.get(bid) or (np.nan, np.nan))
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
        margin=dict(l=60, r=20, t=50, b=60),
        autosize=True, dragmode='pan',
    )
    return fig


def make_source_table(bin_data):
    rows = []
    for ym, entry in sorted(bin_data['months'].items()):
        for s in entry['sources']:
            date_str = str(s['date'])
            rows.append({
                'Year-Month': f"{ym[0]}-{ym[1]:02d}",
                'Date':       f"{date_str[:4]}-{date_str[4:6]}-{date_str[6:8]}",
                'Source':     s.get('source', ''),
                'Qual':       s.get('qual', ''),
                'Lon':        f"{s['lon']:.3f}",
                'Lat':        f"{s['lat']:.3f}",
                'File':       s['file'],
                'Prof index': s['prof_idx'],
            })
    return rows

# ==============================================================================
# App layout
# ==============================================================================

app = Dash(__name__)
app.title = 'Multi-source Anomaly Globe'

DARK = '#1a1a1a'; MID = '#252525'; TEXT = '#dddddd'; ACC = '#88ccff'; DIM = '#aaaaaa'
SCALE_MIN, SCALE_MAX = 1.0, 6.0
VIEW_GLOBE = 'globe'; VIEW_FLAT = 'flat'

_SOURCE_OPTIONS = ([{'label': ' All (combined)', 'value': SOURCE_ALL}]
                   + [{'label': f' {s}', 'value': s} for s in SOURCES])

_DEPTH_MARKS = {i: {'label': f'{depth[i]:.0f} m', 'style': {'color': '#aaa', 'fontSize': '10px'}}
                for i in range(0, N_DEPTH, max(1, N_DEPTH // 10))}

app.index_string = """<!DOCTYPE html>
<html>
<head>
    {%metas%}<title>{%title%}</title>{%favicon%}{%css%}
    <style>
        .zoom-box .rc-slider-rail  { background-color: #555 !important; width: 6px !important; }
        .zoom-box .rc-slider-track { background-color: #88ccff !important; width: 6px !important; }
        .zoom-box .rc-slider-handle {
            border-color: #88ccff !important; background-color: #eee !important;
            width: 16px !important; height: 16px !important; margin-left: -5px !important;
        }
        .zoom-box .rc-slider-handle:hover { border-color: #adf !important; }
    </style>
</head>
<body>{%app_entry%}<footer>{%config%}{%scripts%}{%renderer%}</footer></body>
</html>"""

app.clientside_callback(
    """
    function(_) {
        const LAT_CAP=85, DEG_PER_PX=0.4, CLICK_TOL=5;
        function setup() {
            const globe=document.getElementById('globe');
            if(!globe){setTimeout(setup,200);return;}
            const inner=globe.querySelector('.js-plotly-plot');
            if(!inner||!inner._fullLayout||!window.Plotly){setTimeout(setup,200);return;}
            let dragging=false,moved=false,startX=0,startY=0,startLon=0,startLat=0;
            function currentRotation(){const geo=inner._fullLayout.geo;const r=(geo&&geo.projection&&geo.projection.rotation)||{};return{lon:r.lon||0,lat:r.lat||0};}
            inner.addEventListener('mousedown',function(e){dragging=true;moved=false;startX=e.clientX;startY=e.clientY;const rot=currentRotation();startLon=rot.lon;startLat=rot.lat;e.preventDefault();},true);
            window.addEventListener('mousemove',function(e){if(!dragging)return;const dx=e.clientX-startX,dy=e.clientY-startY;if(Math.abs(dx)+Math.abs(dy)>CLICK_TOL)moved=true;let newLon=startLon-dx*DEG_PER_PX,newLat=startLat+dy*DEG_PER_PX;newLon=((newLon+180)%360+360)%360-180;newLat=Math.max(-LAT_CAP,Math.min(LAT_CAP,newLat));Plotly.relayout(inner,{'geo.projection.rotation.lon':newLon,'geo.projection.rotation.lat':newLat,'geo.projection.rotation.roll':0});},true);
            window.addEventListener('mouseup',function(){dragging=false;},true);
            inner.addEventListener('click',function(e){if(moved){e.stopImmediatePropagation();}},true);
            const SCALE_MIN=1.0,SCALE_MAX=6.0,WHEEL_STEP=0.0015,GESTURE_GAP_MS=400;
            let gz=null,gestureTimer=null,pendingDelta=0,rafId=null;
            function findProjection(){const fl=inner._fullLayout||{};const cands=[fl.geo&&fl.geo._subplot&&fl.geo._subplot.projection,fl._plots&&fl._plots.geo&&fl._plots.geo.projection];const sub=fl.geo&&fl.geo._subplot;if(sub){for(const k in sub){const v=sub[k];if(typeof v==='function'&&typeof v.invert==='function')cands.push(v);}}for(const c of cands){if(c&&typeof c.invert==='function')return c;}return null;}
            const LIMB_CAP_DEG=78,D2R=Math.PI/180;
            function centralAngleDeg(cLat,cLon,lat,lon){const a=cLat*D2R,b=lat*D2R,dl=(lon-cLon)*D2R,c=Math.sin(a)*Math.sin(b)+Math.cos(a)*Math.cos(b)*Math.cos(dl);return Math.acos(Math.max(-1,Math.min(1,c)))/D2R;}
            function cursorLonLat(e){try{const proj=findProjection();if(!proj||!proj.invert)return null;const rect=inner.getBoundingClientRect();const ll=proj.invert([e.clientX-rect.left,e.clientY-rect.top]);if(!ll||!isFinite(ll[0])||!isFinite(ll[1]))return null;const rot=currentRotation();if(centralAngleDeg(rot.lat,rot.lon,ll[1],ll[0])>LIMB_CAP_DEG)return null;return{lon:ll[0],lat:ll[1]};}catch(err){return null;}}
            function wrap180(x){return((x+180)%360+360)%360-180;}
            function flushZoom(){rafId=null;if(pendingDelta===0)return;const geo=inner._fullLayout.geo;const cur=(geo&&geo.projection&&geo.projection.scale)||1.0;let next=cur*Math.exp(-pendingDelta*WHEEL_STEP);next=Math.max(SCALE_MIN,Math.min(SCALE_MAX,next));pendingDelta=0;const relayout={'geo.projection.scale':next};if(gz&&gz.targetLon!==null){const span=SCALE_MAX-gz.startScale;let lin=span>1e-6?(next-gz.startScale)/span:0;lin=Math.max(0,Math.min(1,lin));const _a=0.04,_s=Math.sqrt(_a),f=(Math.sqrt(_a+(1-_a)*lin)-_s)/(1-_s);const dLon=wrap180(gz.targetLon-gz.startLon);relayout['geo.projection.rotation.lon']=wrap180(gz.startLon+f*dLon);relayout['geo.projection.rotation.lat']=gz.startLat+f*(gz.targetLat-gz.startLat);relayout['geo.projection.rotation.roll']=0;}Plotly.relayout(inner,relayout);if(window.dash_clientside&&window.dash_clientside.set_props)window.dash_clientside.set_props('zoom-slider',{value:Math.round(next*10)/10});}
            inner.addEventListener('wheel',function(e){e.preventDefault();e.stopPropagation();if(!gz){const ll=cursorLonLat(e);const tLat=ll?Math.max(-LAT_CAP,Math.min(LAT_CAP,ll.lat)):null;const rot=currentRotation();gz={startScale:(inner._fullLayout.geo.projection.scale)||1.0,startLon:rot.lon,startLat:rot.lat,targetLon:ll?ll.lon:null,targetLat:tLat};}pendingDelta+=e.deltaY;if(rafId===null)rafId=requestAnimationFrame(flushZoom);if(gestureTimer)clearTimeout(gestureTimer);gestureTimer=setTimeout(function(){gz=null;},GESTURE_GAP_MS);},{passive:false,capture:true});
            inner.on('plotly_doubleclick',function(){return false;});
        }
        setup();
        return window.dash_clientside.no_update;
    }
    """,
    Output('globe', 'id'),
    Input('globe', 'id'),
)

app.layout = html.Div(
    style={'backgroundColor': DARK, 'color': TEXT, 'fontFamily': 'sans-serif',
           'height': '100vh', 'display': 'flex', 'flexDirection': 'column',
           'padding': '8px', 'gap': '8px'},
    children=[

    # Header
    html.Div(style={'backgroundColor': MID, 'padding': '8px 14px', 'borderRadius': '4px',
                    'display': 'flex', 'alignItems': 'center', 'gap': '20px'}, children=[
        html.H2('Multi-source Anomaly — Interactive Bin Explorer',
                style={'margin': 0, 'fontSize': '16px', 'color': ACC}),
        html.Span(
            ('Pick a source and color mode, then click a bin.  '
             + (f'Data: {DATASET_FIRST_YM[0]}-{DATASET_FIRST_YM[1]:02d} – '
                f'{DATASET_LAST_YM[0]}-{DATASET_LAST_YM[1]:02d}  |  '
                if DATASET_FIRST_YM and DATASET_LAST_YM else '')
             + f'{len(bins):,} bins  |  {", ".join(SOURCES)}'),
            style={'fontSize': '12px', 'color': DIM}),
    ]),

    # Map row
    html.Div(style={'display': 'flex', 'gap': '8px', 'flex': '0 0 45vh', 'minHeight': '0'},
             children=[

        # Map container
        html.Div(id='map-container',
                 style={'flex': '1', 'backgroundColor': MID, 'borderRadius': '4px',
                        'padding': '4px', 'minHeight': '0', 'position': 'relative'},
                 children=[
            dcc.Graph(id='globe', figure=make_globe(),
                      config={'scrollZoom': False, 'doubleClick': False,
                              'displaylogo': False, 'displayModeBar': False},
                      style={'height': '100%', 'display': 'block'}),
            dcc.Graph(id='flatmap', figure=make_flat(),
                      config={'scrollZoom': True, 'doubleClick': False,
                              'displaylogo': False, 'displayModeBar': False},
                      style={'height': '100%', 'display': 'none'}),
            html.Button('⌂', id='flat-home-btn',
                        style={'display': 'none', 'position': 'absolute',
                               'top': '12px', 'right': '12px', 'zIndex': 10,
                               'fontSize': '18px', 'lineHeight': '1',
                               'padding': '4px 8px', 'cursor': 'pointer',
                               'backgroundColor': 'rgba(20,20,20,0.8)',
                               'color': '#ccc', 'border': '1px solid #666',
                               'borderRadius': '4px'}),
            html.Div(id='zoom-box-wrap', className='zoom-box',
                     style={'position': 'absolute', 'top': '12px', 'left': '12px',
                            'zIndex': 10, 'padding': '10px 8px 14px 8px',
                            'backgroundColor': 'rgba(20,20,20,0.8)',
                            'border': '1px solid #666', 'borderRadius': '6px',
                            'display': 'flex', 'flexDirection': 'column',
                            'alignItems': 'center', 'gap': '6px'},
                     children=[
                html.Div('zoom', style={'fontSize': '10px', 'color': '#ccc', 'fontWeight': 'bold'}),
                dcc.Slider(id='zoom-slider', min=SCALE_MIN, max=SCALE_MAX, step=0.1, value=1.0,
                           vertical=True, verticalHeight=130,
                           marks={SCALE_MIN: {'label': '1×', 'style': {'color': '#ccc'}},
                                  SCALE_MAX: {'label': f'{int(SCALE_MAX)}×', 'style': {'color': '#ccc'}}},
                           tooltip={'placement': 'right', 'always_visible': False}),
            ]),
        ]),

        # Controls panel
        html.Div(style={'width': '230px', 'backgroundColor': MID, 'borderRadius': '4px',
                        'padding': '14px', 'display': 'flex', 'flexDirection': 'column',
                        'gap': '14px', 'overflowY': 'auto', 'minHeight': '0'},
                 children=[

            html.Div([
                html.Label('View', style={'fontSize': '12px', 'color': DIM,
                                          'marginBottom': '6px', 'display': 'block'}),
                dcc.RadioItems(id='view-radio',
                               options=[{'label': ' Globe',    'value': VIEW_GLOBE},
                                        {'label': ' Flat map', 'value': VIEW_FLAT}],
                               value=VIEW_GLOBE,
                               labelStyle={'display': 'inline-block', 'fontSize': '12px',
                                           'marginRight': '12px', 'cursor': 'pointer', 'color': TEXT},
                               inputStyle={'marginRight': '4px'}),
            ]),

            html.Div([
                html.Label('Grid resolution', style={'fontSize': '12px', 'color': DIM,
                                                      'marginBottom': '6px', 'display': 'block'}),
                dcc.RadioItems(id='res-radio',
                               options=[{'label': ' Fine (10242 bins)',  'value': RES_FINE},
                                        {'label': ' Coarse (2562 bins)', 'value': RES_COARSE}],
                               value=RES_FINE,
                               labelStyle={'display': 'block', 'fontSize': '12px',
                                           'marginBottom': '5px', 'cursor': 'pointer', 'color': TEXT},
                               inputStyle={'marginRight': '6px'}),
            ]),

            html.Div([
                html.Label('Source', style={'fontSize': '12px', 'color': DIM,
                                             'marginBottom': '6px', 'display': 'block'}),
                dcc.Dropdown(id='source-select', options=_SOURCE_OPTIONS, value=SOURCE_ALL,
                             clearable=False, style={'fontSize': '12px', 'color': '#111'}),
            ]),

            html.Div([
                html.Label('Color by', style={'fontSize': '12px', 'color': DIM,
                                               'marginBottom': '6px', 'display': 'block'}),
                dcc.Dropdown(id='cmode-select', options=CMODE_OPTIONS, value=CMODE_PROB,
                             clearable=False, style={'fontSize': '12px', 'color': '#111'}),
            ]),

            # Depth slider — shown only when cmode != prob
            html.Div(id='depth-slider-container',
                     style={'display': 'none'},
                     children=[
                html.Label('Depth level', style={'fontSize': '12px', 'color': DIM,
                                                  'marginBottom': '6px', 'display': 'block'}),
                html.Div(id='depth-label',
                         style={'fontSize': '11px', 'color': ACC, 'marginBottom': '4px'}),
                dcc.Slider(id='depth-slider', min=0, max=N_DEPTH - 1, step=1, value=0,
                           marks=_DEPTH_MARKS,
                           tooltip={'placement': 'right', 'always_visible': False}),
            ]),

            html.Div([
                html.Label('Probability threshold', style={'fontSize': '12px', 'color': DIM,
                                                            'marginBottom': '6px', 'display': 'block'}),
                dcc.RadioItems(id='prob-radio',
                               options=[{'label': f' {k}', 'value': k} for k in PROB_KEYS],
                               value='>=1',
                               labelStyle={'display': 'block', 'fontSize': '12px',
                                           'marginBottom': '5px', 'cursor': 'pointer', 'color': TEXT},
                               inputStyle={'marginRight': '6px'}),
            ]),

            html.Div([
                html.Label('Probability basis', style={'fontSize': '12px', 'color': DIM,
                                                        'marginBottom': '6px', 'display': 'block'}),
                dcc.RadioItems(id='basis-radio',
                               options=[{'label': ' Full dataset span',    'value': BASIS_DAY0},
                                        {'label': " Since bin's first obs", 'value': BASIS_SINCE_FIRST}],
                               value=BASIS_DAY0,
                               labelStyle={'display': 'block', 'fontSize': '12px',
                                           'marginBottom': '5px', 'cursor': 'pointer', 'color': TEXT},
                               inputStyle={'marginRight': '6px'}),
            ]),

            html.Div([
                html.Label('Field', style={'fontSize': '12px', 'color': DIM,
                                            'marginBottom': '6px', 'display': 'block'}),
                dcc.RadioItems(id='field-radio',
                               options=[{'label': ' Temperature (T)', 'value': 'T'},
                                        {'label': ' Salinity (S)',    'value': 'S'}],
                               value='T',
                               labelStyle={'display': 'block', 'fontSize': '12px',
                                           'marginBottom': '5px', 'cursor': 'pointer', 'color': TEXT},
                               inputStyle={'marginRight': '6px'}),
            ]),

            html.Div([
                html.Label('Anomaly reference', style={'fontSize': '12px', 'color': DIM,
                                                        'marginBottom': '6px', 'display': 'block'}),
                dcc.RadioItems(id='anom-radio',
                               options=[{'label': ' − climatology',   'value': 'minus_clim'},
                                        {'label': ' − bin time-mean', 'value': 'anom'}],
                               value='minus_clim',
                               labelStyle={'display': 'block', 'fontSize': '12px',
                                           'marginBottom': '5px', 'cursor': 'pointer', 'color': TEXT},
                               inputStyle={'marginRight': '6px'}),
            ]),

            html.Div([
                dcc.Checklist(id='show-sources',
                              options=[{'label': ' Show source table', 'value': 'show'}],
                              value=[],
                              labelStyle={'fontSize': '12px', 'cursor': 'pointer', 'color': TEXT},
                              inputStyle={'marginRight': '6px'}),
            ]),

            html.Div(id='bin-info', style={'fontSize': '11px', 'color': DIM, 'lineHeight': '1.6'}),
        ]),
    ]),

    # Plot row
    html.Div(style={'display': 'flex', 'gap': '8px', 'flex': '1',
                    'minHeight': '0', 'overflow': 'hidden'}, children=[
        html.Div(style={'flex': '1', 'backgroundColor': MID, 'borderRadius': '4px',
                        'padding': '4px', 'minWidth': '0'}, children=[
            dcc.Loading(type='circle', color=ACC, children=[
                dcc.Graph(id='zt-obs', figure=go.Figure(),
                          config={'scrollZoom': True, 'displayModeBar': True,
                                  'modeBarButtonsToAdd': ['resetScale2d'],
                                  'modeBarButtonsToRemove': ['toImage', 'sendDataToCloud'],
                                  'displaylogo': False},
                          style={'height': '100%'}),
            ]),
        ]),
        html.Div(id='zt-anom-container',
                 style={'flex': '1', 'backgroundColor': MID, 'borderRadius': '4px',
                        'padding': '4px', 'minWidth': '0'}, children=[
            dcc.Loading(type='circle', color=ACC, children=[
                dcc.Graph(id='zt-anom', figure=go.Figure(),
                          config={'scrollZoom': True, 'displayModeBar': True,
                                  'modeBarButtonsToAdd': ['resetScale2d'],
                                  'modeBarButtonsToRemove': ['toImage', 'sendDataToCloud'],
                                  'displaylogo': False},
                          style={'height': '100%'}),
            ]),
        ]),
        html.Div(id='source-container',
                 style={'display': 'none', 'flex': '1', 'backgroundColor': MID,
                        'borderRadius': '4px', 'padding': '8px',
                        'flexDirection': 'column', 'overflow': 'hidden', 'minWidth': '0'},
                 children=[
            html.Div('Source profiles', style={'fontSize': '12px', 'color': ACC,
                                               'marginBottom': '6px'}),
            dcc.Loading(type='circle', color=ACC, children=[
                html.Div(id='source-table-container', style={'flex': '1', 'overflowY': 'auto'}),
            ]),
        ]),
    ]),

    dcc.Store(id='selected-bin',  data=None),
    dcc.Store(id='globe-scale',   data=1.0),
    dcc.Store(id='grid-res',      data=RES_FINE),
    dcc.Store(id='map-view',      data=VIEW_GLOBE),
    dcc.Store(id='flat-proj',     data={'scale': 1.0, 'center_lat': 0.0, 'center_lon': 0.0}),
])

# ==============================================================================
# Callbacks
# ==============================================================================

@app.callback(
    Output('depth-slider-container', 'style'),
    Output('depth-label', 'children'),
    Input('cmode-select', 'value'),
    Input('depth-slider', 'value'),
)
def toggle_depth_slider(cmode, depth_idx):
    show = cmode != CMODE_PROB
    style = {'display': 'block' if show else 'none'}
    label = f'{depth[depth_idx or 0]:.1f} m' if show else ''
    return style, label


@app.callback(
    Output('globe', 'style'),
    Output('flatmap', 'style'),
    Output('zoom-box-wrap', 'style'),
    Output('flat-home-btn', 'style'),
    Output('map-view', 'data'),
    Output('flatmap', 'figure', allow_duplicate=True),
    Input('view-radio', 'value'),
    State('source-select', 'value'),
    State('prob-radio',    'value'),
    State('basis-radio',   'value'),
    State('grid-res',      'data'),
    State('selected-bin',  'data'),
    State('cmode-select',  'value'),
    State('depth-slider',  'value'),
    prevent_initial_call=True,
)
def on_view_toggle(view, source_sel, prob_label, basis, res, selected_bin, cmode, depth_idx):
    globe_style = {'height': '100%', 'display': 'block' if view == VIEW_GLOBE else 'none'}
    flat_style  = {'height': '100%', 'display': 'block' if view == VIEW_FLAT  else 'none'}
    zoom_style  = {'position': 'absolute', 'top': '12px', 'left': '12px', 'zIndex': 10,
                   'padding': '10px 8px 14px 8px', 'backgroundColor': 'rgba(20,20,20,0.8)',
                   'border': '1px solid #666', 'borderRadius': '6px',
                   'display': 'flex' if view == VIEW_GLOBE else 'none',
                   'flexDirection': 'column', 'alignItems': 'center', 'gap': '6px'}
    home_style  = {'display': 'block' if view == VIEW_FLAT else 'none',
                   'position': 'absolute', 'top': '12px', 'right': '12px', 'zIndex': 10,
                   'fontSize': '18px', 'lineHeight': '1', 'padding': '4px 8px',
                   'cursor': 'pointer', 'backgroundColor': 'rgba(20,20,20,0.8)',
                   'color': '#ccc', 'border': '1px solid #666', 'borderRadius': '4px'}
    if view == VIEW_FLAT:
        flat_fig = make_flat(source_sel, selected_bin, PROB_KEYS[prob_label], basis,
                             res=res or RES_FINE, cmode=cmode or CMODE_PROB,
                             depth_idx=depth_idx or 0)
    else:
        flat_fig = no_update
    return globe_style, flat_style, zoom_style, home_style, view, flat_fig


@app.callback(
    Output('globe', 'figure', allow_duplicate=True),
    Output('globe-scale', 'data'),
    Input('zoom-slider', 'drag_value'),
    prevent_initial_call=True,
)
def on_zoom(scale):
    scale = scale or 1.0
    patched = Patch()
    patched['layout']['geo']['projection']['scale'] = scale
    return patched, scale


@app.callback(
    Output('globe',        'figure', allow_duplicate=True),
    Output('flatmap',      'figure', allow_duplicate=True),
    Output('grid-res',     'data'),
    Output('selected-bin', 'data',   allow_duplicate=True),
    Input('source-select', 'value'),
    Input('prob-radio',    'value'),
    Input('basis-radio',   'value'),
    Input('res-radio',     'value'),
    Input('cmode-select',  'value'),
    Input('depth-slider',  'value'),
    State('selected-bin',  'data'),
    State('globe-scale',   'data'),
    State('map-view',      'data'),
    prevent_initial_call=True,
)
def on_globe_controls(source_sel, prob_label, basis, res, cmode, depth_idx,
                      selected_bin, scale, view):
    from dash import ctx
    n_threshold = PROB_KEYS[prob_label]
    res         = res or RES_FINE
    cmode       = cmode or CMODE_PROB
    depth_idx   = depth_idx or 0

    if ctx.triggered_id == 'res-radio':
        globe_fig = make_globe(source_sel, None, n_threshold, basis,
                               scale=scale or 1.0, res=res,
                               cmode=cmode, depth_idx=depth_idx)
        flat_fig  = make_flat(source_sel, None, n_threshold, basis, res=res,
                              cmode=cmode, depth_idx=depth_idx)
        return globe_fig, flat_fig, res, None

    a = _globe_arrays(source_sel, n_threshold, basis, cmode, depth_idx, res)
    patched = Patch()
    patched['data'][0]['locations']                 = a['ids']
    patched['data'][0]['z']                         = a['z']
    patched['data'][0]['customdata']                = a['ids']
    patched['data'][0]['colorscale']                = a['colorscale']
    patched['data'][0]['zmin']                      = a['zmin']
    patched['data'][0]['zmax']                      = a['zmax']
    patched['data'][0]['colorbar']['title']['text'] = a['cb_title']
    patched['data'][0]['text']                      = a['hover']
    return patched, patched, res, no_update


@app.callback(
    Output('selected-bin', 'data'),
    Output('globe',   'figure'),
    Output('flatmap', 'figure', allow_duplicate=True),
    Input('globe',   'clickData'),
    Input('flatmap', 'clickData'),
    State('selected-bin', 'data'),
    State('source-select', 'value'),
    State('grid-res', 'data'),
    prevent_initial_call=True,
)
def on_map_click(globe_click, flat_click, current_bin, source_sel, res):
    from dash import ctx
    click_data = globe_click if ctx.triggered_id == 'globe' else flat_click
    if not click_data:
        return current_bin, no_update, no_update
    points = click_data.get('points', [])
    if not points:
        return current_bin, no_update, no_update
    p   = points[0]
    bid = p.get('customdata', p.get('location'))
    b_dict, _, _, _ = _grid(res or RES_FINE)
    if bid is None or bid not in b_dict:
        return current_bin, no_update, no_update
    patched = Patch()
    sel = _selected_trace(bid, res or RES_FINE)
    patched['data'][1]['lon'] = list(sel['lon'])
    patched['data'][1]['lat'] = list(sel['lat'])
    return bid, patched, patched


@app.callback(
    Output('flatmap', 'figure', allow_duplicate=True),
    Output('flat-proj', 'data'),
    Input('flatmap', 'relayoutData'),
    State('flat-proj', 'data'),
    prevent_initial_call=True,
)
def clamp_flat_zoom(relayout, proj_state):
    if not relayout:
        return no_update, no_update
    ps = proj_state or {'scale': 1.0, 'center_lat': 0.0, 'center_lon': 0.0}
    raw_scale = float(relayout.get('geo.projection.scale', ps['scale']))
    raw_lat   = float(relayout.get('geo.center.lat',       ps['center_lat']))
    raw_lon   = float(relayout.get('geo.center.lon',       ps['center_lon']))
    scale      = max(1.0, raw_scale)
    max_clat   = max(0.0, 81.0 - 81.0 / scale)
    center_lat = max(-max_clat, min(max_clat, raw_lat))
    new_state  = {'scale': scale, 'center_lat': center_lat, 'center_lon': raw_lon}
    if scale == raw_scale and center_lat == raw_lat:
        return no_update, new_state
    patched = Patch()
    if scale != raw_scale:
        patched['layout']['geo']['projection']['scale'] = scale
    if center_lat != raw_lat:
        patched['layout']['geo']['center']['lat'] = center_lat
    return patched, new_state


@app.callback(
    Output('flatmap', 'figure', allow_duplicate=True),
    Output('flat-proj', 'data', allow_duplicate=True),
    Input('flat-home-btn', 'n_clicks'),
    State('source-select', 'value'),
    State('prob-radio',    'value'),
    State('basis-radio',   'value'),
    State('grid-res',      'data'),
    State('selected-bin',  'data'),
    State('cmode-select',  'value'),
    State('depth-slider',  'value'),
    prevent_initial_call=True,
)
def flat_home(n_clicks, source_sel, prob_label, basis, res, selected_bin, cmode, depth_idx):
    fig = make_flat(source_sel or SOURCE_ALL, selected_bin,
                    PROB_KEYS.get(prob_label, 1), basis or BASIS_DAY0,
                    res or RES_FINE, uirevision=f'home-{n_clicks}',
                    cmode=cmode or CMODE_PROB, depth_idx=depth_idx or 0)
    return fig, {'scale': 1.0, 'center_lat': 0.0, 'center_lon': 0.0}


@app.callback(
    Output('zt-anom-container', 'style'),
    Output('source-container',  'style'),
    Input('show-sources', 'value'),
)
def toggle_right_panel(show_sources):
    show = 'show' in (show_sources or [])
    anom_style = {'flex': '1', 'backgroundColor': MID, 'borderRadius': '4px',
                  'padding': '4px', 'minWidth': '0',
                  'display': 'none' if show else 'block'}
    src_style  = {'flex': '1', 'backgroundColor': MID, 'borderRadius': '4px',
                  'padding': '8px', 'minWidth': '0', 'overflow': 'hidden',
                  'flexDirection': 'column', 'display': 'flex' if show else 'none'}
    return anom_style, src_style


def _empty_fig(message=''):
    fig = go.Figure()
    fig.update_layout(
        paper_bgcolor='#1a1a1a', plot_bgcolor='#111',
        xaxis=dict(visible=False), yaxis=dict(visible=False),
        annotations=[dict(text=message, xref='paper', yref='paper',
                          x=0.5, y=0.5, showarrow=False,
                          font=dict(size=14, color=DIM))],
        autosize=True, margin=dict(l=20, r=20, t=20, b=20),
    )
    return fig


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


@app.callback(
    Output('zt-anom', 'figure', allow_duplicate=True),
    Input('zt-obs', 'relayoutData'),
    prevent_initial_call=True,
)
def sync_obs_to_anom(relayout):
    p = _axis_patch_from_relayout(relayout)
    return p if p is not None else no_update


@app.callback(
    Output('zt-obs', 'figure', allow_duplicate=True),
    Input('zt-anom', 'relayoutData'),
    prevent_initial_call=True,
)
def sync_anom_to_obs(relayout):
    p = _axis_patch_from_relayout(relayout)
    return p if p is not None else no_update


@app.callback(
    Output('zt-obs',                   'figure'),
    Output('zt-anom',                  'figure'),
    Output('source-table-container',   'children'),
    Output('bin-info',                 'children'),
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
        empty = _empty_fig('Click a bin on the globe')
        return (empty, _empty_fig(''),
                html.Div('No bin selected.', style={'color': DIM, 'fontSize': '12px'}),
                '')

    resolved = resolve_bin(b_dict[bin_id], source_sel, res)
    lon, lat = c_dict.get(bin_id, (np.nan, np.nan))
    n_months = len(resolved['months'])
    n_profs  = resolved['n_profiles']

    if n_months == 0:
        src_note = 'all sources' if source_sel == SOURCE_ALL else source_sel
        msg = _empty_fig(f'Bin {bin_id} ({res_label}) has no {src_note} data')
        return (msg, _empty_fig(''),
                html.Div(f'No {src_note} profiles in this bin.',
                         style={'color': DIM, 'fontSize': '12px'}),
                html.Div(f'Bin {bin_id}: no {src_note} data', style={'color': DIM}))

    obs_fig  = make_zt_figure(resolved, field, depth)
    anom_var = f'{field}_minus_clim' if anom_ref == 'minus_clim' else f'{field}_anom'
    anom_fig = make_zt_figure(resolved, anom_var, depth)

    rows  = make_source_table(resolved)
    table = dash_table.DataTable(
        data=rows,
        columns=[{'name': c, 'id': c} for c in
                 ('Year-Month', 'Date', 'Source', 'Qual', 'Lon', 'Lat', 'File', 'Prof index')],
        style_table={'overflowY': 'auto', 'height': '100%'},
        style_cell={'backgroundColor': '#1e1e1e', 'color': '#ccc', 'fontSize': '11px',
                    'padding': '3px 7px', 'border': '1px solid #333', 'textAlign': 'left',
                    'whiteSpace': 'nowrap', 'overflow': 'hidden',
                    'textOverflow': 'ellipsis', 'maxWidth': '180px'},
        style_header={'backgroundColor': '#2d2d2d', 'color': ACC,
                      'fontWeight': 'bold', 'fontSize': '11px'},
        style_data_conditional=[{'if': {'row_index': 'odd'}, 'backgroundColor': '#232323'}],
        page_size=200, sort_action='native', filter_action='native',
    )

    n_threshold = PROB_KEYS[prob_label]
    prob_val    = bin_probability(b_dict[bin_id], source_sel, n_threshold, basis)
    basis_note  = 'since first obs' if basis == BASIS_SINCE_FIRST else 'full span'
    src_note    = 'All sources' if source_sel == SOURCE_ALL else source_sel
    n_srcs      = len(b_dict[bin_id]['by_source'])
    info = [
        html.Div(f'Bin {bin_id} ({res_label})', style={'color': ACC, 'fontWeight': 'bold'}),
        html.Div(src_note, style={'color': DIM}),
        html.Div(f'lon {lon:.2f}°'),
        html.Div(f'lat {lat:.2f}°'),
        html.Div(f'{n_months} months'),
        html.Div(f'{n_profs} profiles'),
        html.Div(f'{n_srcs} source(s) in bin', style={'color': DIM}),
        html.Div(f'P(≥{n_threshold}) = {prob_val:.3f}'),
        html.Div(f'({basis_note})', style={'color': DIM}),
    ]
    return obs_fig, anom_fig, table, info


if __name__ == '__main__':
    print("Starting app — open http://127.0.0.1:8050")
    app.run(debug=False, host='127.0.0.1', port=8050)
