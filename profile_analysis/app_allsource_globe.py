"""
Interactive multi-source globe app.

Generalizes app_pfl_globe.py to every NCEI-processed source (CTD_WOD, GLD_WOD,
MEOP, MRB_WOD, PFL, PFL_BGC, XBT_WOD, ...) via a Source selector. Reads the
zarr+JSON data store from build_allsource_bin_timeseries.py.

Run with:
    python app_allsource_globe.py
Then open http://127.0.0.1:8050 in a browser.

Controls:
  - Rotate the globe by dragging; zoom with the slider
  - Source selector: "All (combined)" or a single instrument
  - Click a bin to load its Z-T heatmaps and source table
  - Field T / S, anomaly reference (− climatology / − bin time-mean)
  - Probability threshold (>=1/2/3) and denominator basis

Data is PARTITIONED BY SOURCE. The combined ("All") view is derived here on
the fly by profile-weighted pooling of the per-source monthly means (weight =
#profiles that source contributed that month, applied per depth level so a
level missing in one source is filled by the others).

Startup loads only the JSON (metadata + scalars). Per-bin arrays (T/S/Tclim/
Sclim) are read lazily from the zarr store on each bin click.

Requires: dash, plotly, zarr
Data store produced by: build_allsource_bin_timeseries.py
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

_OUTPUT_BASE = Path('/Users/brucel/ecco/yip/profile_data/z_profile_file_analysis/allsource_bin_timeseries')

def _latest_store(base):
    """Return the most recent timestamped store directory, falling back to base itself."""
    parent = base.parent
    candidates = sorted(
        p for p in parent.glob(base.name + '_????????_??????')
        if p.is_dir()
    )
    return candidates[-1] if candidates else base

_app_parser = argparse.ArgumentParser(description='Allsource globe app.')
_app_parser.add_argument('--data', default=None,
                         help='Path to the store directory to load. '
                              'Defaults to the most recently timestamped allsource_bin_timeseries_* directory.')
_app_args, _unknown = _app_parser.parse_known_args()

_DATA_DIR     = Path(_app_args.data) if _app_args.data else _latest_store(_OUTPUT_BASE)
ZARR_STORE    = _DATA_DIR / (_DATA_DIR.name + '.zarr')
METADATA_JSON = _DATA_DIR / (_DATA_DIR.name + '_metadata.json')
print(f"Using data store: {_DATA_DIR}", flush=True)
GEODESIC_2562_CSV = Path('/Users/brucel/ecco/yip/sample_data/ecco-insitu/sweet_gdrive/geodesic/02562_bin_locations.csv')
# Used as fallback if metadata predates the dual-grid rebuild

SOURCE_ALL = 'ALL'   # sentinel for the combined view

# d3-geo (Plotly's geo engine) decides a spherical polygon's interior from ring
# winding, and its convention is the OPPOSITE of the GeoJSON RFC. If cells fill
# their complement (whole globe solid), flip this sign.
#   +1 : ring wound so area-vector points OUTWARD (GeoJSON/CCW-from-outside)
#   -1 : the opposite (what d3-geo usually wants)
CELL_WINDING = -1

# ==============================================================================
# Load data — JSON metadata (fast) + zarr store (lazy arrays on bin click)
# ==============================================================================

print("Loading metadata JSON...", flush=True)
with open(METADATA_JSON, 'r') as _f:
    _meta = json.load(_f)

depth               = np.array(_meta['depth_list'], dtype=np.float32)
SOURCES             = _meta['sources']
DATASET_FIRST_YM    = tuple(_meta['dataset_first_ym']) if _meta.get('dataset_first_ym') else None
DATASET_LAST_YM     = tuple(_meta['dataset_last_ym'])  if _meta.get('dataset_last_ym')  else None
DATASET_MONTHS      = _meta.get('total_months', 1)

# bin_centres: int keys, (lon, lat) tuples
bin_centres         = {int(k): tuple(v) for k, v in _meta['bin_centres'].items()}
bin_centres_coarse  = {int(k): tuple(v) for k, v in _meta['bin_centres_coarse'].items()}

# bins / bins_coarse: lightweight dicts with scalars + profile records only.
# Arrays (T, S, etc.) are NOT loaded here — fetched lazily from zarr on click.
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

# Open zarr store (read-only, lazy — no data is read until indexed)
print("Opening zarr store...", flush=True)
_ZARR = zarr.open_group(str(ZARR_STORE), mode='r')
print("  Ready.")

# Open SQLite profile-records db (read-only; profile rows fetched lazily on bin click)
_PROFILES_DB = _DATA_DIR / (_DATA_DIR.name + '_profiles.db')
_DB_CON = sqlite3.connect(f'file:{_PROFILES_DB}?mode=ro', uri=True, check_same_thread=False)
_DB_CON.row_factory = sqlite3.Row


def _query_profiles(resolution, bid):
    """Return list of profile record dicts for a bin, keyed by (year, month, src)."""
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

N_DEPTH = len(depth)


def _zarr_src_monthly(resolution, bid, src):
    """Return (ym_index, arrays_dict) for one source from the 2D monthly zarr arrays.

    ym_index : (n_months,) int32  — YYYYMM integers, one per row
    arrays_dict : key -> (n_months, N_DEPTH) float32
    """
    node     = _ZARR[f"{resolution}/{bid}/by_source/{src}"]
    ym_index = np.asarray(_ZARR[f"{resolution}/{bid}/ym_index"], dtype=np.int32)
    arrays   = {k: np.asarray(node[k], dtype=np.float32) for k in ('T', 'S', 'Tclim', 'Sclim')}
    return ym_index, arrays


def _zarr_arrays(resolution, bid, ym, src):
    """Read one month's per-source T/S/Tclim/Sclim. Returns dict of (N_DEPTH,) float32."""
    ym_int   = ym[0] * 100 + ym[1]
    ym_index, arrays = _zarr_src_monthly(resolution, bid, src)
    idx = int(np.searchsorted(ym_index, ym_int))
    if idx >= len(ym_index) or ym_index[idx] != ym_int:
        nan = np.full(N_DEPTH, np.nan, dtype=np.float32)
        return {k: nan for k in ('T', 'S', 'Tclim', 'Sclim')}
    return {k: arrays[k][idx] for k in ('T', 'S', 'Tclim', 'Sclim')}


def _zarr_means(resolution, bid, src=None):
    """Read T_mean / S_mean for a bin (combined or per-source) from zarr."""
    path = f"{resolution}/{bid}" if src is None else f"{resolution}/{bid}/by_source/{src}"
    node = _ZARR[path]
    return (np.asarray(node['T_mean'], dtype=np.float32),
            np.asarray(node['S_mean'], dtype=np.float32))

# ==============================================================================
# Geodesic cell polygons (computed once) — the geodesic partitioning as filled
# cells instead of dots. SphericalVoronoi on the bin centres yields the exact
# spherical tiling (12 pentagons + hexagons), data-independent.
# ==============================================================================

def _build_cell_geojson(centres):
    """Return (geojson_dict, ordered_ids). One GeoJSON Polygon feature per bin,
    id = bin_id, vertices in lon/lat. Longitudes are unwrapped relative to each
    cell's centroid so a cell's ring stays contiguous (no ±180 seam splitting the
    polygon); the orthographic projection clips to the visible hemisphere itself.
    """
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
        region = sv.regions[k]
        verts3d = sv.vertices[region]              # (m,3) unit-sphere vertices

        # --- enforce winding so d3/Plotly fills the CELL, not its complement ---
        # On a sphere a ring bounds two areas; orientation picks which is the
        # interior. Reverse the ring unless it is counter-clockwise as seen from
        # OUTSIDE at the cell centre, i.e. its area-vector points outward.
        area_vec = np.zeros(3)
        for i in range(len(verts3d)):
            area_vec += np.cross(verts3d[i], verts3d[(i + 1) % len(verts3d)])
        # CELL_WINDING selects which orientation d3-geo treats as "interior".
        if np.dot(area_vec, xyz[k]) * CELL_WINDING < 0:   # xyz[k] = outward normal
            region = region[::-1]

        ring = [xyz2lonlat(sv.vertices[i]) for i in region]
        lons = np.array([p[0] for p in ring])
        lats = np.array([p[1] for p in ring])
        # Unwrap longitudes to be continuous around this cell's own centre lon,
        # so a cell straddling ±180 doesn't wrap the wrong way across the globe.
        c_lon = centres[bid][0]
        lons = c_lon + ((lons - c_lon + 180.0) % 360.0 - 180.0)
        coords = [[float(a), float(b)] for a, b in zip(lons, lats)]
        coords.append(coords[0])   # close the ring
        features.append({
            'type': 'Feature',
            'id': int(bid),
            'geometry': {'type': 'Polygon', 'coordinates': [coords]},
            'properties': {},
        })
    return {'type': 'FeatureCollection', 'features': features}, ids


print("Building fine (10242-bin) geodesic cell polygons...", flush=True)
CELL_GEOJSON, CELL_IDS = _build_cell_geojson(bin_centres)
print(f"  {len(CELL_IDS):,} fine cells built")

print("Building coarse (2562-bin) geodesic cell polygons...", flush=True)
# Always build from the full 2562-centre CSV so every cell has a polygon,
# even if only a subset of bins have data in the store.
_coarse_raw = np.genfromtxt(GEODESIC_2562_CSV, delimiter=',')
_all_coarse_centres = {i + 1: (float(_coarse_raw[i, 0]), float(_coarse_raw[i, 1]))
                       for i in range(len(_coarse_raw))}
# Merge into bin_centres_coarse so hover/click lookups cover all cells
if not bin_centres_coarse:
    bin_centres_coarse = _all_coarse_centres
else:
    for k, v in _all_coarse_centres.items():
        bin_centres_coarse.setdefault(k, v)
COARSE_GEOJSON, COARSE_IDS = _build_cell_geojson(_all_coarse_centres)
print(f"  {len(COARSE_IDS):,} coarse cells built")


def _months_between(ym0, ym1):
    """Inclusive count of calendar months from ym0 to ym1, each (year, month)."""
    return (ym1[0] - ym0[0]) * 12 + (ym1[1] - ym0[1]) + 1

# Threshold selector maps label -> N
PROB_KEYS = {'>=1': 1, '>=2': 2, '>=3': 3}

# Denominator-basis selector
BASIS_DAY0        = 'day0'          # full dataset span (same denom for all bins)
BASIS_SINCE_FIRST = 'since_first'   # bin's first obs -> dataset end

# ==============================================================================
# Per-source resolution
# ==============================================================================

def _pool_level_weighted(entries, key):
    """Profile-weighted pool of a per-source monthly mean across `entries`.

    entries: list of by_source dicts, each with key -> (N_DEPTH,) and 'n'.
    Weighting is applied per depth level, so a level that is NaN for one source
    is still filled from the sources that have it. For a single-element list
    this returns that source's own mean unchanged.
    """
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


RES_FINE   = 'fine'    # 10242-bin grid
RES_COARSE = 'coarse'  # 2562-bin grid


def _source_entries(month, source_sel):
    """List of by_source entries for the selected source (or all) in a month."""
    bs = month['by_source']
    if source_sel == SOURCE_ALL:
        return list(bs.values())
    entry = bs.get(source_sel)
    return [entry] if entry is not None else []


def resolve_bin(bin_data, source_sel, res=RES_FINE):
    """Collapse a partitioned bin to a flat shape for the given source (or All).

    Returns a dict shaped like the PFL-era bin format so make_zt_figure and
    make_source_table work unchanged:
      {'bin_id', 'months': {ym: {'T','S','Tclim','Sclim','sources','n'}},
       'T_mean','S_mean','hit_counts','first_ym','last_ym', 'n_profiles'}

    Arrays are fetched lazily from zarr here (one zarr read per bin click, not
    at startup).
    """
    bid  = bin_data['bin_id']
    zres = 'fine' if res == RES_FINE else 'coarse'

    # Fetch profile records for this bin from SQLite (one query for all months/sources)
    db_recs = _query_profiles(zres, bid)

    months = {}
    for ym, m in bin_data['months'].items():
        if source_sel != SOURCE_ALL and source_sel not in m['by_source']:
            continue
        src_list = (list(m['by_source'].keys()) if source_sel == SOURCE_ALL
                    else [source_sel])
        entries = []
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
        hit_counts = bin_data['hit_counts']
        first_ym   = bin_data['first_ym']
        last_ym    = bin_data['last_ym']
        n_profiles = sum(r['n_profiles'] for r in bin_data['by_source'].values())
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


def bin_has_source(bin_data, source_sel):
    """True if this bin has any data for the selected source (or All)."""
    return source_sel == SOURCE_ALL or source_sel in bin_data['by_source']


def bin_probability(bin_data, source_sel, n_threshold, basis):
    """Probability that the bin received >= n_threshold profiles in a month for
    the selected source (or combined), under the chosen denominator basis."""
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

# ==============================================================================
# Globe figure
# ==============================================================================

def _grid(res):
    """Return (bins_dict, centres_dict, geojson, ids) for the chosen resolution."""
    if res == RES_COARSE:
        return bins_coarse, bin_centres_coarse, COARSE_GEOJSON, COARSE_IDS
    return bins, bin_centres, CELL_GEOJSON, CELL_IDS


def _visible_ids(source_sel, res=RES_FINE):
    b_dict, c_dict, _, _ = _grid(res)
    return [b for b in sorted(b_dict.keys())
            if b in c_dict and bin_has_source(b_dict[b], source_sel)]


def _build_hover(ids_clean, source_sel, n_threshold, basis, prob_label, res=RES_FINE):
    b_dict, c_dict, _, _ = _grid(res)
    src_note = 'All sources' if source_sel == SOURCE_ALL else source_sel
    return [
        f'Bin {b} — {src_note}'
        f'<br>lon={c_dict[b][0]:.2f}, lat={c_dict[b][1]:.2f}'
        f'<br>P({prob_label} prof/month) = {bin_probability(b_dict[b], source_sel, n_threshold, basis):.3f}'
        for b in ids_clean
    ]


def _globe_arrays(source_sel, n_threshold, basis, res=RES_FINE):
    b_dict, _, _, _ = _grid(res)
    ids_clean  = _visible_ids(source_sel, res)
    probs      = [bin_probability(b_dict[b], source_sel, n_threshold, basis) for b in ids_clean]
    prob_label = next(k for k, v in PROB_KEYS.items() if v == n_threshold)
    hover      = _build_hover(ids_clean, source_sel, n_threshold, basis, prob_label, res)
    return dict(ids=ids_clean, probs=probs, hover=hover, prob_label=prob_label)


def _selected_trace(selected_bin_id, res=RES_FINE):
    """A thin Scattergeo outline of the selected cell."""
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


def make_globe(source_sel=SOURCE_ALL, selected_bin_id=None,
               n_threshold=1, basis=BASIS_DAY0, scale=1.0, res=RES_FINE):
    _, _, geojson, _ = _grid(res)
    a = _globe_arrays(source_sel, n_threshold, basis, res)

    fig = go.Figure(go.Choropleth(
        geojson=geojson,
        locations=a['ids'],
        z=a['probs'],
        featureidkey='id',
        colorscale='RdYlGn',
        zmin=0, zmax=1,
        marker=dict(line=dict(color='rgba(40,40,40,0.5)', width=0.3)),
        colorbar=dict(title=f'P({a["prob_label"]})', thickness=14, len=0.6),
        text=a['hover'],
        hoverinfo='text',
        customdata=a['ids'],
    ))
    fig.add_trace(_selected_trace(selected_bin_id, res))

    fig.update_geos(
        projection_type='orthographic',
        projection_scale=scale,
        projection_rotation=dict(lon=0, lat=20, roll=0),
        showland=True,    landcolor='#2a2a2a',
        showocean=True,   oceancolor='#0d1b2a',
        showlakes=True,   lakecolor='#3a3a3a',
        showrivers=True,  rivercolor='#3a3a3a',
        showcoastlines=True, coastlinecolor='#cccccc',
        coastlinewidth=0.8,
        showframe=False,
        bgcolor='#1a1a1a',
    )
    fig.update_layout(
        margin=dict(l=0, r=0, t=0, b=0),
        paper_bgcolor='#1a1a1a',
        geo=dict(bgcolor='#1a1a1a'),
        autosize=True,
        dragmode=False,
        uirevision='globe',
    )
    return fig


def _build_coastline_trace():
    """Scattergeo trace of 50m coastlines — rendered on top of the choropleth."""
    lons, lats = [], []
    for g in cfeature.NaturalEarthFeature('physical', 'coastline', '50m').geometries():
        parts = [g] if g.geom_type == 'LineString' else list(g.geoms)
        for part in parts:
            x, y = part.xy
            lons.extend(list(x) + [None])
            lats.extend(list(y) + [None])
    return go.Scattergeo(
        lon=lons, lat=lats,
        mode='lines',
        line=dict(color='#4a9eff', width=0.9),
        hoverinfo='skip',
        showlegend=False,
    )

_COAST_TRACE = _build_coastline_trace()


def make_flat(source_sel=SOURCE_ALL, selected_bin_id=None,
              n_threshold=1, basis=BASIS_DAY0, res=RES_FINE, uirevision='flat'):
    """Flat (natural-earth projection) choropleth — same data as make_globe."""
    _, _, geojson, _ = _grid(res)
    a = _globe_arrays(source_sel, n_threshold, basis, res)

    fig = go.Figure(go.Choropleth(
        geojson=geojson,
        locations=a['ids'],
        z=a['probs'],
        featureidkey='id',
        colorscale='RdYlGn',
        zmin=0, zmax=1,
        marker=dict(line=dict(color='rgba(40,40,40,0.5)', width=0.3)),
        colorbar=dict(title=f'P({a["prob_label"]})', thickness=14, len=0.6),
        text=a['hover'],
        hoverinfo='text',
        customdata=a['ids'],
    ))
    fig.add_trace(_selected_trace(selected_bin_id, res))  # trace[1]
    fig.add_trace(_COAST_TRACE)                            # trace[2]

    fig.update_geos(
        projection_type='natural earth',
        showland=True,    landcolor='#2a2a2a',
        showocean=True,   oceancolor='#0d1b2a',
        showlakes=True,   lakecolor='#3a3a3a',
        showrivers=True,  rivercolor='#3a3a3a',
        showcoastlines=False,
        showframe=False,
        bgcolor='#1a1a1a',
        lataxis_range=[-90, 90],
        lonaxis_range=[-180, 180],
    )
    fig.update_layout(
        margin=dict(l=0, r=0, t=0, b=0),
        paper_bgcolor='#1a1a1a',
        geo=dict(bgcolor='#1a1a1a'),
        autosize=True,
        dragmode='pan',
        uirevision=uirevision,
    )
    return fig

# ==============================================================================
# Z-T heatmap helpers  (operate on a RESOLVED bin — flat month shape)
# ==============================================================================

def make_zt_figure(bin_data, var, depth):
    """
    var: 'T' | 'S' | 'Tclim' | 'Sclim' | 'T_minus_clim' | 'S_minus_clim'
         | 'T_anom' | 'S_anom'
    """
    months_dict = bin_data['months']
    sorted_yms  = sorted(months_dict.keys())
    if not sorted_yms:
        return _empty_fig('No data for this source in this bin')

    x = [f'{ym[0]}-{ym[1]:02d}-15' for ym in sorted_yms]

    is_diff_clim = var in ('T_minus_clim', 'S_minus_clim')
    is_anom      = var in ('T_anom', 'S_anom')
    base_letter  = var[0]   # 'T' or 'S'

    if is_diff_clim:
        obs_key  = base_letter
        clim_key = base_letter + 'clim'
        matrix = np.stack(
            [months_dict[ym][obs_key] - months_dict[ym][clim_key] for ym in sorted_yms],
            axis=0).T
        colorscale = 'RdBu_r'
        abs_max = np.nanpercentile(np.abs(matrix[~np.isnan(matrix)]), 97) if not np.all(np.isnan(matrix)) else 1
        zmin, zmax = -abs_max, abs_max
        title_var = f'{base_letter} − climatology'
        units = '°C' if base_letter == 'T' else 'PSU'

    elif is_anom:
        matrix = np.stack([months_dict[ym][base_letter] for ym in sorted_yms], axis=0).T
        long_mean = bin_data[f'{base_letter}_mean']
        matrix = matrix - long_mean[:, np.newaxis]
        colorscale = 'RdBu_r'
        abs_max = np.nanpercentile(np.abs(matrix[~np.isnan(matrix)]), 97) if not np.all(np.isnan(matrix)) else 1
        zmin, zmax = -abs_max, abs_max
        title_var = f'{base_letter} − bin time-mean'
        units = '°C' if base_letter == 'T' else 'PSU'

    else:
        matrix = np.stack([months_dict[ym][var] for ym in sorted_yms], axis=0).T
        valid  = matrix[~np.isnan(matrix)]
        if base_letter == 'T':
            colorscale = 'Thermal'
            zmin = np.nanpercentile(valid, 2)  if len(valid) else -2
            zmax = np.nanpercentile(valid, 98) if len(valid) else 30
        else:
            colorscale = 'Haline'
            zmin = np.nanpercentile(valid, 2)  if len(valid) else 30
            zmax = np.nanpercentile(valid, 98) if len(valid) else 38
        clim_suffix = ' (clim)' if var.endswith('clim') else ''
        title_var = base_letter + clim_suffix
        units = '°C' if base_letter == 'T' else 'PSU'

    fig = go.Figure(go.Heatmap(
        x=x, y=depth, z=matrix,
        colorscale=colorscale, zmin=zmin, zmax=zmax,
        colorbar=dict(title=units, thickness=12),
        hovertemplate='%{x|%b %Y}<br>Depth: %{y} m<br>Value: %{z:.3f}<extra></extra>',
    ))

    bid = bin_data['bin_id']
    lon, lat = (bin_centres.get(bid) or bin_centres_coarse.get(bid) or (np.nan, np.nan))
    fig.update_layout(
        title=dict(
            text=f'Bin {bid} ({lon:.1f}°, {lat:.1f}°) — {title_var} (monthly bin average)',
            font=dict(size=13, color='#ddd'),
        ),
        xaxis=dict(title='', color='#aaa', gridcolor='#333', automargin=True,
                   type='date', tickformat='%Y', tickformatstops=[
                       dict(dtickrange=[None, 'M3'],  value='%b %Y'),
                       dict(dtickrange=['M3', 'M18'], value='%b %Y'),
                       dict(dtickrange=['M18', None], value='%Y'),
                   ]),
        yaxis=dict(title='Depth (m)', autorange='reversed', color='#aaa', gridcolor='#333'),
        paper_bgcolor='#1a1a1a',
        plot_bgcolor='#111',
        margin=dict(l=60, r=20, t=50, b=60),
        autosize=True,
        dragmode='pan',   # left-drag pans; wheel still zooms (scrollZoom=True)
    )
    return fig


def make_source_table(bin_data):
    rows = []
    for ym, entry in sorted(bin_data['months'].items()):
        for s in entry['sources']:
            date_str = str(s['date'])
            date_fmt = f"{date_str[:4]}-{date_str[4:6]}-{date_str[6:8]}"
            rows.append({
                'Year-Month': f"{ym[0]}-{ym[1]:02d}",
                'Date':       date_fmt,
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
app.title = 'Multi-source Profile Globe'

app.index_string = """<!DOCTYPE html>
<html>
<head>
    {%metas%}<title>{%title%}</title>{%favicon%}{%css%}
    <style>
        .zoom-box .rc-slider-rail  { background-color: #555 !important; width: 6px !important; }
        .zoom-box .rc-slider-track { background-color: #88ccff !important; width: 6px !important; }
        .zoom-box .rc-slider-handle {
            border-color: #88ccff !important;
            background-color: #eee !important;
            width: 16px !important; height: 16px !important;
            margin-left: -5px !important;
        }
        .zoom-box .rc-slider-handle:hover { border-color: #adf !important; }
    </style>
</head>
<body>
    {%app_entry%}
    <footer>{%config%}{%scripts%}{%renderer%}</footer>
</body>
</html>"""

# Custom globe rotation (see app_pfl_globe.py for rationale).
app.clientside_callback(
    """
    function(_) {
        const LAT_CAP    = 85;
        const DEG_PER_PX = 0.4;
        const CLICK_TOL  = 5;

        function setup() {
            const globe = document.getElementById('globe');
            if (!globe) { setTimeout(setup, 200); return; }
            const inner = globe.querySelector('.js-plotly-plot');
            if (!inner || !inner._fullLayout || !window.Plotly) {
                setTimeout(setup, 200); return;
            }

            let dragging = false, moved = false;
            let startX = 0, startY = 0, startLon = 0, startLat = 0;

            function currentRotation() {
                const geo = inner._fullLayout.geo;
                const r = (geo && geo.projection && geo.projection.rotation) || {};
                return {lon: r.lon || 0, lat: r.lat || 0};
            }

            inner.addEventListener('mousedown', function(e) {
                dragging = true; moved = false;
                startX = e.clientX; startY = e.clientY;
                const rot = currentRotation();
                startLon = rot.lon; startLat = rot.lat;
                e.preventDefault();
            }, true);

            window.addEventListener('mousemove', function(e) {
                if (!dragging) return;
                const dx = e.clientX - startX;
                const dy = e.clientY - startY;
                if (Math.abs(dx) + Math.abs(dy) > CLICK_TOL) moved = true;

                let newLon = startLon - dx * DEG_PER_PX;
                let newLat = startLat + dy * DEG_PER_PX;
                newLon = ((newLon + 180) % 360 + 360) % 360 - 180;
                newLat = Math.max(-LAT_CAP, Math.min(LAT_CAP, newLat));

                Plotly.relayout(inner, {
                    'geo.projection.rotation.lon': newLon,
                    'geo.projection.rotation.lat': newLat,
                    'geo.projection.rotation.roll': 0,
                });
            }, true);

            window.addEventListener('mouseup', function() { dragging = false; }, true);

            inner.addEventListener('click', function(e) {
                if (moved) { e.stopImmediatePropagation(); }
            }, true);

            // Wheel zoom-toward-cursor. We adjust ONLY scale + rotation via
            // Plotly.relayout (never roll), deliberately bypassing Plotly's
            // built-in geo scrollZoom, which couples lon/lat/roll/scale into one
            // gesture and glitches over the poles (the original bug).
            //
            // Behaviour: on the first wheel tick of a gesture we capture the
            // lon/lat under the cursor and the starting rotation/scale. As scale
            // climbs from start toward SCALE_MAX we ease the rotation from its
            // start value toward "cursor-centred", so there is no jump at gesture
            // start and the point is fully centred at max zoom. Zooming back out
            // un-eases it. Falls back to plain centre-locked zoom if the
            // projection's invert() isn't reachable.
            const SCALE_MIN = 1.0, SCALE_MAX = 6.0;
            const WHEEL_STEP = 0.0015;    // per wheel delta unit; feels natural
            const GESTURE_GAP_MS = 400;   // idle gap that ends a zoom gesture

            let gz = null;   // active gesture state, or null between gestures
            let gestureTimer = null;

            // Locate the d3 geo projection Plotly builds for the subplot. The
            // handle is internal/undocumented and has moved between versions, so
            // probe several known locations and remember the first that inverts.
            function findProjection() {
                const fl = inner._fullLayout || {};
                const cands = [
                    fl.geo && fl.geo._subplot && fl.geo._subplot.projection,
                    fl.geo && fl.geo._subplot && fl.geo._subplot.mockAxis,  // placeholder
                    fl._plots && fl._plots.geo && fl._plots.geo.projection,
                    inner._fullData && fl.geo && fl.geo.projection,
                ];
                // Deep fallback: scan fullLayout.geo._subplot for any fn with .invert
                const sub = fl.geo && fl.geo._subplot;
                if (sub) {
                    for (const k in sub) {
                        const v = sub[k];
                        if (typeof v === 'function' && typeof v.invert === 'function') {
                            cands.push(v);
                        }
                    }
                }
                for (const c of cands) {
                    if (c && typeof c.invert === 'function') return c;
                }
                return null;
            }

            // Reject cursor targets within this angular distance of the limb.
            // Orthographic's inverse Jacobian blows up at the rim (d angle / d
            // pixel -> infinity), so a cursor a few px from the edge inverts to a
            // wildly imprecise, jittery point. Inside ~78deg of the view centre
            // the mapping is well-conditioned; beyond it we fall back to plain
            // centre-locked zoom instead of chasing a bad point.
            const LIMB_CAP_DEG = 78;
            const D2R = Math.PI / 180;

            function centralAngleDeg(centerLat, centerLon, lat, lon) {
                const a = centerLat * D2R, b = lat * D2R;
                const dl = (lon - centerLon) * D2R;
                const c = Math.sin(a) * Math.sin(b)
                        + Math.cos(a) * Math.cos(b) * Math.cos(dl);
                return Math.acos(Math.max(-1, Math.min(1, c))) / D2R;
            }

            function cursorLonLat(e) {
                try {
                    const proj = findProjection();
                    if (!proj || !proj.invert) return null;
                    const rect = inner.getBoundingClientRect();
                    const px = e.clientX - rect.left;
                    const py = e.clientY - rect.top;
                    const ll = proj.invert([px, py]);
                    if (!ll || !isFinite(ll[0]) || !isFinite(ll[1])) return null;
                    // Reject near-limb points (unstable invert). Plotly's
                    // rotation.lon is the centre longitude in degrees East, so
                    // the current view centre is (lat = rot.lat, lon = rot.lon).
                    const rot = currentRotation();
                    const ang = centralAngleDeg(rot.lat, rot.lon, ll[1], ll[0]);
                    if (ang > LIMB_CAP_DEG) return null;
                    return {lon: ll[0], lat: ll[1]};
                } catch (err) { return null; }
            }

            function wrap180(x) { return ((x + 180) % 360 + 360) % 360 - 180; }

            // Wheel deltas are coalesced and applied ONCE per animation frame via
            // requestAnimationFrame. Trackpads/mice fire wheel events far faster
            // than Plotly can re-render; relayout-per-event queues up and stutters.
            // Accumulating deltaY and flushing at ~60fps smooths the motion.
            let pendingDelta = 0;
            let rafId = null;

            function flushZoom() {
                rafId = null;
                if (pendingDelta === 0) return;
                const geo = inner._fullLayout.geo;
                const cur = (geo && geo.projection && geo.projection.scale) || 1.0;

                // multiplicative so each notch is a constant zoom ratio
                let next = cur * Math.exp(-pendingDelta * WHEEL_STEP);
                next = Math.max(SCALE_MIN, Math.min(SCALE_MAX, next));
                pendingDelta = 0;

                const relayout = {'geo.projection.scale': next};

                if (gz && gz.targetLon !== null) {
                    // Fraction of the way from gesture-start scale to max zoom,
                    // then a concave ease (sqrt, normalised to hit 0 at lin=0 and
                    // 1 at lin=1, with finite slope at 0 so the first tiny scroll
                    // doesn't kick). Front-loaded so the pan is visible early.
                    const span = SCALE_MAX - gz.startScale;
                    let lin = span > 1e-6 ? (next - gz.startScale) / span : 0;
                    lin = Math.max(0, Math.min(1, lin));
                    const _a = 0.04, _s = Math.sqrt(_a);
                    const f = (Math.sqrt(_a + (1 - _a) * lin) - _s) / (1 - _s);
                    // Shortest-path longitude interpolation; latitude is linear.
                    const dLon = wrap180(gz.targetLon - gz.startLon);
                    relayout['geo.projection.rotation.lon'] = wrap180(gz.startLon + f * dLon);
                    relayout['geo.projection.rotation.lat'] = gz.startLat + f * (gz.targetLat - gz.startLat);
                    relayout['geo.projection.rotation.roll'] = 0;
                }

                Plotly.relayout(inner, relayout);

                // Sync the slider handle without triggering on_zoom: that callback
                // listens to `drag_value`, so setting `value` moves the handle but
                // does NOT loop back into a relayout.
                if (window.dash_clientside && window.dash_clientside.set_props) {
                    window.dash_clientside.set_props('zoom-slider',
                        {value: Math.round(next * 10) / 10});
                }
            }

            inner.addEventListener('wheel', function(e) {
                e.preventDefault();          // stop the page from scrolling
                e.stopPropagation();

                // Start a new gesture if none active (first tick after a pause).
                if (!gz) {
                    const ll = cursorLonLat(e);
                    // Clamp target latitude to the same +/-LAT_CAP the drag
                    // handler respects, so centring on a near-pole point can't
                    // drive rotation into the unstable polar zone (flip/wobble).
                    const tLat = ll ? Math.max(-LAT_CAP, Math.min(LAT_CAP, ll.lat)) : null;
                    const rot = currentRotation();
                    gz = {
                        startScale: (inner._fullLayout.geo.projection.scale) || 1.0,
                        startLon:   rot.lon,
                        startLat:   rot.lat,
                        // rotation.lon is centre-longitude in degrees East, so to
                        // bring a point of longitude L to centre we set lon = L.
                        targetLon:  ll ? ll.lon : null,
                        targetLat:  tLat,
                    };
                }

                // Accumulate this event's delta; flush on the next frame.
                // (Scale + slider sync both happen in flushZoom.)
                pendingDelta += e.deltaY;
                if (rafId === null) rafId = requestAnimationFrame(flushZoom);

                // Reset the gesture after a short idle gap so the next scroll
                // re-captures a fresh cursor target.
                if (gestureTimer) clearTimeout(gestureTimer);
                gestureTimer = setTimeout(function() { gz = null; }, GESTURE_GAP_MS);
            }, {passive: false, capture: true});

            inner.on('plotly_doubleclick', function() { return false; });
        }
        setup();
        return window.dash_clientside.no_update;
    }
    """,
    Output('globe', 'id'),
    Input('globe', 'id'),
)

DARK = '#1a1a1a'
MID  = '#252525'
TEXT = '#dddddd'
ACC  = '#88ccff'
DIM  = '#aaaaaa'   # secondary / label text — readable on dark background

SCALE_MIN, SCALE_MAX = 1.0, 6.0

VIEW_GLOBE = 'globe'
VIEW_FLAT  = 'flat'

_SOURCE_OPTIONS = ([{'label': ' All (combined)', 'value': SOURCE_ALL}]
                   + [{'label': f' {s}', 'value': s} for s in SOURCES])

app.layout = html.Div(style={'backgroundColor': DARK, 'color': TEXT,
                              'fontFamily': 'sans-serif', 'height': '100vh',
                              'display': 'flex', 'flexDirection': 'column',
                              'padding': '8px', 'gap': '8px'}, children=[

    html.Div(style={'backgroundColor': MID, 'padding': '8px 14px',
                    'borderRadius': '4px', 'display': 'flex',
                    'alignItems': 'center', 'gap': '20px'}, children=[
        html.H2('Multi-source — Interactive Bin Explorer',
                style={'margin': 0, 'fontSize': '16px', 'color': ACC}),
        html.Span(
            ('Pick a source, then click a bin to inspect its Z-T structure.  '
             + (f'Data: {DATASET_FIRST_YM[0]}-{DATASET_FIRST_YM[1]:02d} – '
                f'{DATASET_LAST_YM[0]}-{DATASET_LAST_YM[1]:02d}  |  '
                if DATASET_FIRST_YM and DATASET_LAST_YM else '')
             + f'{len(bins):,} bins  |  {", ".join(SOURCES)}'),
            style={'fontSize': '12px', 'color': DIM}),
    ]),

    html.Div(style={'display': 'flex', 'gap': '8px', 'flex': '0 0 45vh',
                    'minHeight': '0'}, children=[

        html.Div(id='map-container',
                 style={'flex': '1', 'backgroundColor': MID,
                        'borderRadius': '4px', 'padding': '4px',
                        'minHeight': '0', 'position': 'relative'}, children=[
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
                            'alignItems': 'center', 'gap': '6px'}, children=[
                html.Div('zoom', style={'fontSize': '10px', 'color': '#ccc',
                                        'fontWeight': 'bold'}),
                dcc.Slider(
                    id='zoom-slider',
                    min=SCALE_MIN, max=SCALE_MAX, step=0.1, value=1.0,
                    vertical=True, verticalHeight=130,
                    marks={SCALE_MIN: {'label': '1×', 'style': {'color': '#ccc'}},
                           SCALE_MAX: {'label': f'{int(SCALE_MAX)}×', 'style': {'color': '#ccc'}}},
                    tooltip={'placement': 'right', 'always_visible': False},
                ),
            ]),
        ]),

        html.Div(style={'width': '210px', 'backgroundColor': MID,
                        'borderRadius': '4px', 'padding': '14px',
                        'display': 'flex', 'flexDirection': 'column',
                        'gap': '14px', 'overflowY': 'auto', 'minHeight': '0'},
                 children=[

            html.Div([
                html.Label('View',
                           style={'fontSize': '12px', 'color': DIM,
                                  'marginBottom': '6px', 'display': 'block'}),
                dcc.RadioItems(
                    id='view-radio',
                    options=[
                        {'label': ' Globe',    'value': VIEW_GLOBE},
                        {'label': ' Flat map', 'value': VIEW_FLAT},
                    ],
                    value=VIEW_FLAT,
                    labelStyle={'display': 'inline-block', 'fontSize': '12px',
                                'marginRight': '12px', 'cursor': 'pointer',
                                'color': TEXT},
                    inputStyle={'marginRight': '4px'},
                ),
            ]),

            html.Div([
                html.Label('Grid resolution',
                           style={'fontSize': '12px', 'color': DIM,
                                  'marginBottom': '6px', 'display': 'block'}),
                dcc.RadioItems(
                    id='res-radio',
                    options=[
                        {'label': ' Fine (10242 bins)',   'value': RES_FINE},
                        {'label': ' Coarse (2562 bins)',  'value': RES_COARSE},
                    ],
                    value=RES_FINE,
                    labelStyle={'display': 'block', 'fontSize': '12px',
                                'marginBottom': '5px', 'cursor': 'pointer',
                                'color': TEXT},
                    inputStyle={'marginRight': '6px'},
                ),
            ]),

            html.Div([
                html.Label('Source',
                           style={'fontSize': '12px', 'color': DIM,
                                  'marginBottom': '6px', 'display': 'block'}),
                dcc.Dropdown(
                    id='source-select',
                    options=_SOURCE_OPTIONS,
                    value=SOURCE_ALL,
                    clearable=False,
                    style={'fontSize': '12px', 'color': '#111'},
                ),
            ]),

            html.Div([
                html.Label('Globe: probability threshold',
                           style={'fontSize': '12px', 'color': DIM,
                                  'marginBottom': '6px', 'display': 'block'}),
                dcc.RadioItems(
                    id='prob-radio',
                    options=[{'label': f' {k}', 'value': k} for k in PROB_KEYS],
                    value='>=1',
                    labelStyle={'display': 'block', 'fontSize': '12px',
                                'marginBottom': '5px', 'cursor': 'pointer',
                                'color': TEXT},
                    inputStyle={'marginRight': '6px'},
                ),
            ]),

            html.Div([
                html.Label('Probability basis',
                           style={'fontSize': '12px', 'color': DIM,
                                  'marginBottom': '6px', 'display': 'block'}),
                dcc.RadioItems(
                    id='basis-radio',
                    options=[
                        {'label': ' Full dataset span', 'value': BASIS_DAY0},
                        {'label': ' Since bin\'s first obs', 'value': BASIS_SINCE_FIRST},
                    ],
                    value=BASIS_DAY0,
                    labelStyle={'display': 'block', 'fontSize': '12px',
                                'marginBottom': '5px', 'cursor': 'pointer',
                                'color': TEXT},
                    inputStyle={'marginRight': '6px'},
                ),
            ]),

            html.Div([
                html.Label('Field', style={'fontSize': '12px', 'color': DIM,
                                           'marginBottom': '6px', 'display': 'block'}),
                dcc.RadioItems(
                    id='field-radio',
                    options=[
                        {'label': ' Temperature (T)', 'value': 'T'},
                        {'label': ' Salinity (S)',    'value': 'S'},
                    ],
                    value='T',
                    labelStyle={'display': 'block', 'fontSize': '12px',
                                'marginBottom': '5px', 'cursor': 'pointer',
                                'color': TEXT},
                    inputStyle={'marginRight': '6px'},
                ),
            ]),

            html.Div([
                html.Label('Anomaly reference', style={'fontSize': '12px', 'color': DIM,
                                                       'marginBottom': '6px', 'display': 'block'}),
                dcc.RadioItems(
                    id='anom-radio',
                    options=[
                        {'label': ' − climatology',    'value': 'minus_clim'},
                        {'label': ' − bin time-mean',  'value': 'anom'},
                    ],
                    value='minus_clim',
                    labelStyle={'display': 'block', 'fontSize': '12px',
                                'marginBottom': '5px', 'cursor': 'pointer',
                                'color': TEXT},
                    inputStyle={'marginRight': '6px'},
                ),
            ]),

            html.Div([
                dcc.Checklist(
                    id='show-sources',
                    options=[{'label': ' Show source table', 'value': 'show'}],
                    value=[],
                    labelStyle={'fontSize': '12px', 'cursor': 'pointer', 'color': TEXT},
                    inputStyle={'marginRight': '6px'},
                ),
            ]),

            html.Div(id='bin-info', style={'fontSize': '11px', 'color': DIM,
                                            'lineHeight': '1.6'}),
        ]),
    ]),

    html.Div(style={'display': 'flex', 'gap': '8px', 'flex': '1',
                    'minHeight': '0', 'overflow': 'hidden'}, children=[

        html.Div(style={'flex': '1', 'backgroundColor': MID,
                        'borderRadius': '4px', 'padding': '4px',
                        'minWidth': '0'}, children=[
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
                 style={'flex': '1', 'backgroundColor': MID,
                        'borderRadius': '4px', 'padding': '4px',
                        'minWidth': '0'}, children=[
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
                        'flexDirection': 'column', 'overflow': 'hidden',
                        'minWidth': '0'}, children=[
            html.Div('Source profiles', style={'fontSize': '12px', 'color': ACC,
                                                'marginBottom': '6px'}),
            dcc.Loading(type='circle', color=ACC, children=[
                html.Div(id='source-table-container', style={'flex': '1', 'overflowY': 'auto'}),
            ]),
        ]),
    ]),

    dcc.Store(id='selected-bin', data=None),
    dcc.Store(id='globe-scale', data=1.0),
    dcc.Store(id='grid-res', data=RES_FINE),
    dcc.Store(id='map-view', data=VIEW_FLAT),
    dcc.Store(id='flat-proj', data={'scale': 1.0, 'center_lat': 0.0, 'center_lon': 0.0}),
])

# ==============================================================================
# Callbacks
# ==============================================================================

@app.callback(
    Output('globe',         'style'),
    Output('flatmap',       'style'),
    Output('zoom-box-wrap', 'style'),
    Output('flat-home-btn', 'style'),
    Output('map-view',      'data'),
    Output('flatmap',       'figure', allow_duplicate=True),
    Input('view-radio', 'value'),
    State('source-select', 'value'),
    State('prob-radio',    'value'),
    State('basis-radio',   'value'),
    State('grid-res',      'data'),
    State('selected-bin',  'data'),
    prevent_initial_call=True,
)
def on_view_toggle(view, source_sel, prob_label, basis, res, selected_bin):
    globe_style = {'height': '100%', 'display': 'block' if view == VIEW_GLOBE else 'none'}
    flat_style  = {'height': '100%', 'display': 'block' if view == VIEW_FLAT  else 'none'}
    zoom_style  = {'position': 'absolute', 'top': '12px', 'left': '12px',
                   'zIndex': 10, 'padding': '10px 8px 14px 8px',
                   'backgroundColor': 'rgba(20,20,20,0.8)',
                   'border': '1px solid #666', 'borderRadius': '6px',
                   'display': 'flex' if view == VIEW_GLOBE else 'none',
                   'flexDirection': 'column', 'alignItems': 'center', 'gap': '6px'}
    home_style  = {'display': 'block' if view == VIEW_FLAT else 'none',
                   'position': 'absolute', 'top': '12px', 'right': '12px',
                   'zIndex': 10, 'fontSize': '18px', 'lineHeight': '1',
                   'padding': '4px 8px', 'cursor': 'pointer',
                   'backgroundColor': 'rgba(20,20,20,0.8)',
                   'color': '#ccc', 'border': '1px solid #666', 'borderRadius': '4px'}
    if view == VIEW_FLAT:
        n_threshold = PROB_KEYS[prob_label]
        flat_fig = make_flat(source_sel, selected_bin, n_threshold, basis,
                             res=res or RES_FINE)
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
    Output('globe',   'figure', allow_duplicate=True),
    Output('flatmap', 'figure', allow_duplicate=True),
    Output('grid-res', 'data'),
    Output('selected-bin', 'data', allow_duplicate=True),
    Input('source-select', 'value'),
    Input('prob-radio', 'value'),
    Input('basis-radio', 'value'),
    Input('res-radio', 'value'),
    State('selected-bin', 'data'),
    State('globe-scale', 'data'),
    State('map-view', 'data'),
    prevent_initial_call=True,
)
def on_globe_controls(source_sel, prob_label, basis, res, selected_bin, scale, view):
    from dash import ctx
    n_threshold = PROB_KEYS[prob_label]
    res = res or RES_FINE

    if ctx.triggered_id == 'res-radio':
        globe_fig = make_globe(source_sel, None, n_threshold, basis,
                               scale=scale or 1.0, res=res)
        flat_fig  = make_flat(source_sel, None, n_threshold, basis, res=res)
        return globe_fig, flat_fig, res, None

    a = _globe_arrays(source_sel, n_threshold, basis, res)
    patched = Patch()
    patched['data'][0]['locations']                 = a['ids']
    patched['data'][0]['z']                         = a['probs']
    patched['data'][0]['customdata']                = a['ids']
    patched['data'][0]['colorbar']['title']['text'] = f'P({prob_label})'
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
    if click_data is None:
        return current_bin, no_update, no_update
    points = click_data.get('points', [])
    if not points:
        return current_bin, no_update, no_update
    p = points[0]
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
    center_lon = raw_lon

    new_state = {'scale': scale, 'center_lat': center_lat, 'center_lon': center_lon}

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
    prevent_initial_call=True,
)
def flat_home(n_clicks, source_sel, n_threshold, basis, res, selected_bin):
    # Return a full fresh figure rather than a patch so the entire geo layout
    # (center, projection scale, rotation, axis ranges) resets cleanly.
    # Pass n_clicks as uirevision so Plotly sees a new value every press and
    # does not preserve the prior pan/zoom interaction state.
    fig = make_flat(source_sel or SOURCE_ALL,
                    selected_bin,
                    n_threshold or 1,
                    basis or BASIS_DAY0,
                    res or RES_FINE,
                    uirevision=f'home-{n_clicks}')
    return fig, {'scale': 1.0, 'center_lat': 0.0, 'center_lon': 0.0}


def _empty_fig(message):
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
    patch = Patch()
    touched = False
    if 'xaxis.range[0]' in relayout and 'xaxis.range[1]' in relayout:
        patch['layout']['xaxis']['range'] = [relayout['xaxis.range[0]'],
                                             relayout['xaxis.range[1]']]
        patch['layout']['xaxis']['autorange'] = False
        touched = True
    if 'yaxis.range[0]' in relayout and 'yaxis.range[1]' in relayout:
        patch['layout']['yaxis']['range'] = [relayout['yaxis.range[0]'],
                                             relayout['yaxis.range[1]']]
        patch['layout']['yaxis']['autorange'] = False
        touched = True
    if relayout.get('xaxis.autorange'):
        patch['layout']['xaxis']['autorange'] = True
        touched = True
    if relayout.get('yaxis.autorange'):
        patch['layout']['yaxis']['autorange'] = 'reversed'
        touched = True
    return patch if touched else None


@app.callback(
    Output('zt-anom', 'figure', allow_duplicate=True),
    Input('zt-obs', 'relayoutData'),
    prevent_initial_call=True,
)
def sync_obs_to_anom(relayout):
    patch = _axis_patch_from_relayout(relayout)
    return patch if patch is not None else no_update


@app.callback(
    Output('zt-obs', 'figure', allow_duplicate=True),
    Input('zt-anom', 'relayoutData'),
    prevent_initial_call=True,
)
def sync_anom_to_obs(relayout):
    patch = _axis_patch_from_relayout(relayout)
    return patch if patch is not None else no_update


@app.callback(
    Output('zt-anom-container', 'style'),
    Output('source-container', 'style'),
    Input('show-sources', 'value'),
)
def toggle_right_panel(show_sources):
    show = 'show' in (show_sources or [])
    anom_style = {'flex': '1', 'backgroundColor': MID, 'borderRadius': '4px',
                  'padding': '4px', 'minWidth': '0',
                  'display': 'none' if show else 'block'}
    src_style  = {'flex': '1', 'backgroundColor': MID, 'borderRadius': '4px',
                  'padding': '8px', 'minWidth': '0',
                  'overflow': 'hidden', 'flexDirection': 'column',
                  'display': 'flex' if show else 'none'}
    return anom_style, src_style


@app.callback(
    Output('zt-obs', 'figure'),
    Output('zt-anom', 'figure'),
    Output('source-table-container', 'children'),
    Output('bin-info', 'children'),
    Input('selected-bin', 'data'),
    Input('source-select', 'value'),
    Input('field-radio', 'value'),
    Input('anom-radio', 'value'),
    Input('prob-radio', 'value'),
    Input('basis-radio', 'value'),
    Input('grid-res', 'data'),
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

    rows = make_source_table(resolved)
    table = dash_table.DataTable(
        data=rows,
        columns=[{'name': c, 'id': c} for c in
                 ('Year-Month', 'Date', 'Source', 'Qual', 'Lon', 'Lat', 'File', 'Prof index')],
        style_table={'overflowY': 'auto', 'height': '100%'},
        style_cell={'backgroundColor': '#1e1e1e', 'color': '#ccc',
                    'fontSize': '11px', 'padding': '3px 7px',
                    'border': '1px solid #333', 'textAlign': 'left',
                    'whiteSpace': 'nowrap', 'overflow': 'hidden',
                    'textOverflow': 'ellipsis', 'maxWidth': '180px'},
        style_header={'backgroundColor': '#2d2d2d', 'color': ACC,
                      'fontWeight': 'bold', 'fontSize': '11px'},
        style_data_conditional=[
            {'if': {'row_index': 'odd'}, 'backgroundColor': '#232323'},
        ],
        page_size=200,
        sort_action='native',
        filter_action='native',
    )

    n_threshold = PROB_KEYS[prob_label]
    prob_glyph  = prob_label.replace('>=', '≥')
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
        html.Div(f'P({prob_glyph}) = {prob_val:.3f}'),
        html.Div(f'({basis_note})', style={'color': DIM}),
    ]

    return obs_fig, anom_fig, table, info


if __name__ == '__main__':
    print("Starting app — open http://127.0.0.1:8050")
    app.run(debug=False, host='127.0.0.1', port=8050)
