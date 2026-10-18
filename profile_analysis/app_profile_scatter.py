"""
Profile scatter app — shows individual profile locations on globe/flat map.

Colour-codes each dot by source instrument. A dual range-slider filters the
date window. Clicking any dot loads the enclosing bin's full T/S timeseries
in the lower panels (same zarr-backed data as app_allsource_anomaly_globe.py).

Spatial decimation keeps at most ~20 000 points on screen: cell size is derived
from the actual bounding box of the query result, so a tight cluster uses small
cells (many points kept) while a globally-spread query uses large cells (one
representative per large area).

Requires the data store built by build_allsource_bin_timeseries_with_anomalies.py
(zarr + _profiles.db).

Run with:
    python app_profile_scatter.py
Then open http://127.0.0.1:8050
"""

import json
import math
import sqlite3
import argparse
import numpy as np
import plotly.graph_objects as go
import cartopy.feature as cfeature
import zarr
from dash import Dash, dcc, html, Input, Output, State, no_update, Patch
from pathlib import Path

# ==============================================================================
# Config
# ==============================================================================

_OUTPUT_BASE = Path('/Users/brucel/ecco/yip/profile_data/z_profile_file_analysis/'
                    'allsource_bin_timeseries_with_anomalies')

def _latest_store(base):
    parent = base.parent
    candidates = sorted(p for p in parent.glob(base.name + '_????????_??????') if p.is_dir())
    return candidates[-1] if candidates else base

_app_parser = argparse.ArgumentParser(description='Profile scatter app.')
_app_parser.add_argument('--data', default=None,
                         help='Path to anomaly store directory (default: most recent timestamped '
                              'allsource_bin_timeseries_with_anomalies_* directory).')
_app_args, _unknown = _app_parser.parse_known_args()

_DATA_DIR     = Path(_app_args.data) if _app_args.data else _latest_store(_OUTPUT_BASE)
ZARR_STORE    = _DATA_DIR / (_DATA_DIR.name + '.zarr')
METADATA_JSON = _DATA_DIR / (_DATA_DIR.name + '_metadata.json')
PROFILES_DB   = _DATA_DIR / (_DATA_DIR.name + '_profiles.db')
print(f"Using data store: {_DATA_DIR}", flush=True)

SOURCE_ALL     = 'ALL'
RES_FINE       = 'fine'
RES_COARSE     = 'coarse'
VIEW_GLOBE     = 'globe'
VIEW_FLAT      = 'flat'
TARGET_N       = 20_000
SCALE_MIN, SCALE_MAX = 1.0, 6.0

# ==============================================================================
# Metadata
# ==============================================================================

print("Loading metadata JSON...", flush=True)
with open(METADATA_JSON) as _f:
    _meta = json.load(_f)

depth              = np.array(_meta['depth_list'], dtype=np.float32)
SOURCES            = _meta['sources']
DATASET_FIRST_YM   = tuple(_meta['dataset_first_ym']) if _meta.get('dataset_first_ym') else (1990, 1)
DATASET_LAST_YM    = tuple(_meta['dataset_last_ym'])  if _meta.get('dataset_last_ym')  else (2024, 12)
DATASET_MONTHS     = _meta.get('total_months', 1)
bin_centres        = {int(k): tuple(v) for k, v in _meta['bin_centres'].items()}
bin_centres_coarse = {int(k): tuple(v) for k, v in _meta['bin_centres_coarse'].items()}
N_DEPTH            = len(depth)


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
print(f"  {len(bins):,} fine bins, {len(bins_coarse):,} coarse bins, {len(SOURCES)} sources")

print("Opening zarr store...", flush=True)
_ZARR = zarr.open_group(str(ZARR_STORE), mode='r')
print("  Ready.")

# Read-only SQLite; opened once, queries are thread-safe via check_same_thread=False
_DB_CON = sqlite3.connect(f'file:{PROFILES_DB}?mode=ro', uri=True, check_same_thread=False)
_DB_CON.row_factory = sqlite3.Row

# ==============================================================================
# Source colours (10-colour qualitative palette)
# ==============================================================================

_PALETTE = ['#636EFA', '#EF553B', '#00CC96', '#AB63FA', '#FFA15A',
            '#19D3F3', '#FF6692', '#B6E880', '#FF97FF', '#FECB52']
_SRC_COLOR = {src: _PALETTE[i % len(_PALETTE)] for i, src in enumerate(SOURCES)}

# ==============================================================================
# Spatial decimation
# ==============================================================================

def _decimate(lats, lons, target_n=TARGET_N):
    """Boolean mask keeping ≤target_n spatially spread points.

    Cell size is derived from the bounding box of the data, so a tight
    regional query keeps more points than a global one with the same M.
    """
    n = len(lats)
    if n <= target_n:
        return np.ones(n, dtype=bool)

    lat_range = float(lats.max() - lats.min())
    lon_range = float(lons.max() - lons.min())
    mean_lat  = float(np.mean(lats))

    lat_km   = lat_range * 111.0
    lon_km   = lon_range * 111.0 * max(0.1, math.cos(math.radians(mean_lat)))
    area_km2 = max(lat_km * lon_km, 1.0)

    cell_km  = math.sqrt(area_km2 / target_n)
    cell_deg = max(0.01, cell_km / 111.0)

    lb = np.floor(lats / cell_deg).astype(np.int64)
    bb = np.floor(lons / cell_deg).astype(np.int64)

    # Shuffle so the kept point per cell is random, not always the first in db order
    rng  = np.random.default_rng(seed=0)
    perm = rng.permutation(n)
    keys = lb[perm] * 1_000_000 + bb[perm]
    _, first_in_perm = np.unique(keys, return_index=True)
    selected = perm[first_in_perm]

    mask = np.zeros(n, dtype=bool)
    mask[selected] = True
    return mask


# ==============================================================================
# DB query
# ==============================================================================

def query_profiles(source, res, year_min, year_max):
    """Query and spatially decimate profile locations.

    Returns dict with arrays + counts, or None if no rows found.
    Note: queries by year*100+month; the profiles table has no composite index
    on year/month so this is a full scan filtered post-index on (resolution).
    """
    ym_min = year_min * 100 + 1
    ym_max = year_max * 100 + 12

    if source == SOURCE_ALL:
        cur = _DB_CON.execute(
            'SELECT lat, lon, bin_id, src, year, month FROM profiles '
            'WHERE resolution=? AND (year*100+month) BETWEEN ? AND ?',
            (res, ym_min, ym_max),
        )
    else:
        cur = _DB_CON.execute(
            'SELECT lat, lon, bin_id, src, year, month FROM profiles '
            'WHERE resolution=? AND src=? AND (year*100+month) BETWEEN ? AND ?',
            (res, source, ym_min, ym_max),
        )
    rows = cur.fetchall()
    if not rows:
        return None

    lats    = np.array([r['lat']    for r in rows], dtype=np.float32)
    lons    = np.array([r['lon']    for r in rows], dtype=np.float32)
    bin_ids = np.array([r['bin_id'] for r in rows], dtype=np.int32)
    srcs    = np.array([r['src']    for r in rows])
    dates   = np.array([r['year'] * 100 + r['month'] for r in rows], dtype=np.int32)

    n_total = len(lats)
    mask    = _decimate(lats, lons)

    return {
        'lats':    lats[mask],
        'lons':    lons[mask],
        'bin_ids': bin_ids[mask],
        'srcs':    srcs[mask],
        'dates':   dates[mask],
        'n_shown': int(mask.sum()),
        'n_total': n_total,
    }


# ==============================================================================
# Zarr + bin helpers (mirrors anomaly app)
# ==============================================================================

def _zarr_src_monthly(resolution, bid, src):
    node     = _ZARR[f"{resolution}/{bid}/by_source/{src}"]
    ym_index = np.asarray(_ZARR[f"{resolution}/{bid}/ym_index"], dtype=np.int32)
    arrays   = {k: np.asarray(node[k], dtype=np.float32) for k in ('T', 'S', 'Tclim', 'Sclim')}
    return ym_index, arrays


def _zarr_arrays(resolution, bid, ym, src):
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


def _query_profiles_for_bin(resolution, bid):
    cur = _DB_CON.execute(
        'SELECT src, year, month, date, lon, lat, file, prof_idx, qual '
        'FROM profiles WHERE resolution=? AND bin_id=?',
        (resolution, bid),
    )
    result = {}
    for row in cur:
        ym  = (row['year'], row['month'])
        src = row['src']
        result.setdefault(ym, {}).setdefault(src, []).append({
            'date': row['date'], 'lon': row['lon'], 'lat': row['lat'],
            'file': row['file'], 'prof_idx': row['prof_idx'],
            'qual': row['qual'], 'source': src,
        })
    return result


def _grid(res):
    return (bins_coarse, bin_centres_coarse) if res == RES_COARSE else (bins, bin_centres)


def resolve_bin(bin_data, source_sel, res):
    bid     = bin_data['bin_id']
    db_recs = _query_profiles_for_bin(res, bid)
    months  = {}
    for ym, m in bin_data['months'].items():
        if source_sel != SOURCE_ALL and source_sel not in m['by_source']:
            continue
        src_list = list(m['by_source'].keys()) if source_sel == SOURCE_ALL else [source_sel]
        entries  = []
        for src in src_list:
            if src not in m['by_source']:
                continue
            arrays = _zarr_arrays(res, bid, ym, src)
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
        T_mean, S_mean = _zarr_means(res, bid)
        n_profiles = sum(r['n_profiles'] for r in bin_data['by_source'].values())
    else:
        T_mean, S_mean = _zarr_means(res, bid, source_sel)
        n_profiles = bin_data['by_source'][source_sel]['n_profiles']
    return {
        'bin_id':     bid,
        'months':     months,
        'T_mean':     np.asarray(T_mean, dtype=float),
        'S_mean':     np.asarray(S_mean, dtype=float),
        'n_profiles': n_profiles,
    }


# ==============================================================================
# Map figures
# ==============================================================================

GEO_STYLE = dict(
    showland=True,    landcolor='#2a2a2a',
    showocean=True,   oceancolor='#0d1b2a',
    showlakes=True,   lakecolor='#3a3a3a',
    showrivers=True,  rivercolor='#3a3a3a',
    showcoastlines=True, coastlinecolor='#cccccc', coastlinewidth=0.8,
    showframe=False,  bgcolor='#1a1a1a',
)


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


print("Building coastline trace...", flush=True)
_COAST_TRACE = _build_coastline_trace()


def _scatter_traces(data, source_sel):
    """One Scattergeo trace per source (for legend); empty trace if no data."""
    if data is None:
        return [go.Scattergeo(lat=[], lon=[], mode='markers',
                              showlegend=False, name='_empty')]

    if source_sel == SOURCE_ALL:
        traces = []
        for src in SOURCES:
            mask = data['srcs'] == src
            if not np.any(mask):
                continue
            lats_s    = data['lats'][mask]
            lons_s    = data['lons'][mask]
            bin_ids_s = data['bin_ids'][mask]
            dates_s   = data['dates'][mask]
            hover = [
                f'{src}<br>{d // 100}-{d % 100:02d}'
                f'<br>lon={lo:.2f}, lat={la:.2f}'
                for d, lo, la in zip(dates_s, lons_s, lats_s)
            ]
            traces.append(go.Scattergeo(
                lat=lats_s.tolist(), lon=lons_s.tolist(),
                mode='markers',
                marker=dict(size=3, color=_SRC_COLOR[src], opacity=0.7),
                name=src, text=hover, hoverinfo='text',
                customdata=bin_ids_s.tolist(),
            ))
        return traces or [go.Scattergeo(lat=[], lon=[], mode='markers',
                                        showlegend=False, name='_empty')]
    else:
        hover = [
            f'{s}<br>{d // 100}-{d % 100:02d}'
            f'<br>lon={lo:.2f}, lat={la:.2f}'
            for s, d, lo, la in zip(data['srcs'], data['dates'],
                                    data['lons'], data['lats'])
        ]
        return [go.Scattergeo(
            lat=data['lats'].tolist(), lon=data['lons'].tolist(),
            mode='markers',
            marker=dict(size=3, color=_SRC_COLOR.get(source_sel, '#88ccff'), opacity=0.7),
            name=source_sel, text=hover, hoverinfo='text',
            customdata=data['bin_ids'].tolist(),
        )]


def _selection_trace(sel_lon, sel_lat):
    if sel_lon is None:
        return go.Scattergeo(lat=[], lon=[], mode='markers',
                             hoverinfo='skip', showlegend=False, name='_sel')
    return go.Scattergeo(
        lat=[sel_lat], lon=[sel_lon], mode='markers',
        marker=dict(size=14, symbol='star', color='white',
                    line=dict(color='black', width=1)),
        hoverinfo='skip', showlegend=False, name='_sel',
    )


_LEGEND_STYLE = dict(
    bgcolor='rgba(20,20,20,0.75)', font=dict(color='#ccc', size=11),
    x=0.01, y=0.01, xanchor='left', yanchor='bottom',
)


def make_globe(data, source_sel, sel_lon=None, sel_lat=None, scale=1.0):
    traces = _scatter_traces(data, source_sel) + [_selection_trace(sel_lon, sel_lat)]
    fig = go.Figure(traces)
    fig.update_geos(projection_type='orthographic',
                    projection_scale=scale,
                    projection_rotation=dict(lon=0, lat=20, roll=0),
                    **GEO_STYLE)
    fig.update_layout(margin=dict(l=0, r=0, t=0, b=0),
                      paper_bgcolor='#1a1a1a', geo=dict(bgcolor='#1a1a1a'),
                      autosize=True, dragmode=False, uirevision='globe',
                      legend=_LEGEND_STYLE)
    return fig


def make_flat(data, source_sel, sel_lon=None, sel_lat=None, uirevision='flat'):
    traces = _scatter_traces(data, source_sel) + [_selection_trace(sel_lon, sel_lat),
                                                   _COAST_TRACE]
    fig = go.Figure(traces)
    fig.update_geos(projection_type='natural earth',
                    showcoastlines=False,
                    lataxis_range=[-90, 90], lonaxis_range=[-180, 180],
                    **GEO_STYLE)
    fig.update_layout(margin=dict(l=0, r=0, t=0, b=0),
                      paper_bgcolor='#1a1a1a', geo=dict(bgcolor='#1a1a1a'),
                      autosize=True, dragmode='pan', uirevision=uirevision,
                      legend=_LEGEND_STYLE)
    return fig


# ==============================================================================
# Z-T heatmap (identical logic to anomaly app)
# ==============================================================================

def make_zt_figure(bin_data, var, depth_arr):
    months_dict = bin_data['months']
    sorted_yms  = sorted(months_dict.keys())
    if not sorted_yms:
        return _empty_fig('No data for this source in this bin')
    x            = [f'{ym[0]}-{ym[1]:02d}-15' for ym in sorted_yms]
    is_diff_clim = var in ('T_minus_clim', 'S_minus_clim')
    is_anom      = var in ('T_anom', 'S_anom')
    base_letter  = var[0]
    if is_diff_clim:
        obs_key, clim_key = base_letter, base_letter + 'clim'
        matrix = np.stack([months_dict[ym][obs_key] - months_dict[ym][clim_key]
                           for ym in sorted_yms], axis=0).T
        colorscale = 'RdBu_r'
        abs_max    = (np.nanpercentile(np.abs(matrix[np.isfinite(matrix)]), 97)
                      if np.any(np.isfinite(matrix)) else 1)
        zmin, zmax = -abs_max, abs_max
        title_var, units = f'{base_letter} − climatology', '°C' if base_letter == 'T' else 'PSU'
    elif is_anom:
        matrix = np.stack([months_dict[ym][base_letter] for ym in sorted_yms], axis=0).T
        matrix = matrix - bin_data[f'{base_letter}_mean'][:, np.newaxis]
        colorscale = 'RdBu_r'
        abs_max    = (np.nanpercentile(np.abs(matrix[np.isfinite(matrix)]), 97)
                      if np.any(np.isfinite(matrix)) else 1)
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
        x=x, y=depth_arr, z=matrix,
        colorscale=colorscale, zmin=zmin, zmax=zmax,
        colorbar=dict(title=units, thickness=12),
        hovertemplate='%{x|%b %Y}<br>Depth: %{y} m<br>Value: %{z:.3f}<extra></extra>',
    ))
    bid      = bin_data['bin_id']
    lon, lat = bin_centres.get(bid) or bin_centres_coarse.get(bid) or (np.nan, np.nan)
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


def _empty_fig(message=''):
    fig = go.Figure()
    fig.update_layout(
        paper_bgcolor='#1a1a1a', plot_bgcolor='#111',
        xaxis=dict(visible=False), yaxis=dict(visible=False),
        annotations=[dict(text=message, xref='paper', yref='paper',
                          x=0.5, y=0.5, showarrow=False,
                          font=dict(size=14, color='#aaaaaa'))],
        autosize=True, margin=dict(l=20, r=20, t=20, b=20),
    )
    return fig


# ==============================================================================
# App setup
# ==============================================================================

app = Dash(__name__)
app.title = 'Profile Scatter'

DARK = '#1a1a1a'; MID = '#252525'; TEXT = '#dddddd'; ACC = '#88ccff'; DIM = '#aaaaaa'

_SOURCE_OPTIONS = ([{'label': ' All sources', 'value': SOURCE_ALL}]
                   + [{'label': f' {s}', 'value': s} for s in SOURCES])

YEAR_MIN = DATASET_FIRST_YM[0]
YEAR_MAX = DATASET_LAST_YM[0]
_YEAR_MARKS = {
    y: {'label': str(y), 'style': {'color': '#aaa', 'fontSize': '10px'}}
    for y in range(YEAR_MIN, YEAR_MAX + 1, 5)
}

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

# Globe drag + zoom — identical to anomaly app
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

# ==============================================================================
# Layout
# ==============================================================================

app.layout = html.Div(
    style={'backgroundColor': DARK, 'color': TEXT, 'fontFamily': 'sans-serif',
           'height': '100vh', 'display': 'flex', 'flexDirection': 'column',
           'padding': '8px', 'gap': '8px'},
    children=[

    # Header
    html.Div(style={'backgroundColor': MID, 'padding': '8px 14px', 'borderRadius': '4px',
                    'display': 'flex', 'alignItems': 'center', 'gap': '20px'}, children=[
        html.H2('Profile Scatter — Location & Provenance',
                style={'margin': 0, 'fontSize': '16px', 'color': ACC}),
        html.Span(
            ('Click a dot to see the bin timeseries.  '
             + (f'Data: {DATASET_FIRST_YM[0]}-{DATASET_FIRST_YM[1]:02d} – '
                f'{DATASET_LAST_YM[0]}-{DATASET_LAST_YM[1]:02d}  |  '
                if DATASET_FIRST_YM and DATASET_LAST_YM else '')
             + f'{", ".join(SOURCES)}'),
            style={'fontSize': '12px', 'color': DIM}),
        html.Span(id='count-label', style={'fontSize': '12px', 'color': ACC,
                                            'marginLeft': 'auto'}),
    ]),

    # Map row
    html.Div(style={'display': 'flex', 'gap': '8px', 'flex': '0 0 48vh', 'minHeight': '0'},
             children=[

        # Map container
        html.Div(id='map-container',
                 style={'flex': '1', 'backgroundColor': MID, 'borderRadius': '4px',
                        'padding': '4px', 'minHeight': '0', 'position': 'relative'},
                 children=[
            dcc.Loading(type='circle', color=ACC, children=[
                dcc.Graph(id='globe', figure=make_globe(None, SOURCE_ALL),
                          config={'scrollZoom': False, 'doubleClick': False,
                                  'displaylogo': False, 'displayModeBar': False},
                          style={'height': '100%', 'display': 'none'}),
                dcc.Graph(id='flatmap', figure=make_flat(None, SOURCE_ALL),
                          config={'scrollZoom': True, 'doubleClick': False,
                                  'displaylogo': False, 'displayModeBar': False},
                          style={'height': '100%', 'display': 'block'}),
            ]),
            html.Button('⌂', id='flat-home-btn',
                        style={'display': 'block', 'position': 'absolute',
                               'top': '12px', 'right': '12px', 'zIndex': 10,
                               'fontSize': '18px', 'lineHeight': '1',
                               'padding': '4px 8px', 'cursor': 'pointer',
                               'backgroundColor': 'rgba(20,20,20,0.8)',
                               'color': '#ccc', 'border': '1px solid #666',
                               'borderRadius': '4px'}),
            html.Div(id='zoom-box-wrap', className='zoom-box',
                     style={'display': 'none', 'position': 'absolute',
                            'top': '12px', 'left': '12px', 'zIndex': 10,
                            'padding': '10px 8px 14px 8px',
                            'backgroundColor': 'rgba(20,20,20,0.8)',
                            'border': '1px solid #666', 'borderRadius': '6px',
                            'flexDirection': 'column', 'alignItems': 'center', 'gap': '6px'},
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
                               value=VIEW_FLAT,
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
                html.Label('Time window', style={'fontSize': '12px', 'color': DIM,
                                                  'marginBottom': '6px', 'display': 'block'}),
                html.Div(id='year-range-label',
                         style={'fontSize': '11px', 'color': ACC, 'marginBottom': '4px'}),
                dcc.RangeSlider(
                    id='year-slider',
                    min=YEAR_MIN, max=YEAR_MAX, step=1,
                    value=[YEAR_MIN, YEAR_MAX],
                    marks=_YEAR_MARKS,
                    tooltip={'placement': 'bottom', 'always_visible': False},
                    allowCross=False,
                ),
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

            html.Div(id='bin-info', style={'fontSize': '11px', 'color': DIM, 'lineHeight': '1.6'}),
        ]),
    ]),

    # Lower panels
    html.Div(style={'display': 'flex', 'gap': '8px', 'flex': '1',
                    'minHeight': '0', 'overflow': 'hidden'}, children=[
        html.Div(style={'flex': '1', 'backgroundColor': MID, 'borderRadius': '4px',
                        'padding': '4px', 'minWidth': '0'}, children=[
            dcc.Loading(type='circle', color=ACC, children=[
                dcc.Graph(id='zt-obs', figure=_empty_fig('Click a profile dot on the map'),
                          config={'scrollZoom': True, 'displayModeBar': True,
                                  'modeBarButtonsToAdd': ['resetScale2d'],
                                  'modeBarButtonsToRemove': ['toImage', 'sendDataToCloud'],
                                  'displaylogo': False},
                          style={'height': '100%'}),
            ]),
        ]),
        html.Div(style={'flex': '1', 'backgroundColor': MID, 'borderRadius': '4px',
                        'padding': '4px', 'minWidth': '0'}, children=[
            dcc.Loading(type='circle', color=ACC, children=[
                dcc.Graph(id='zt-anom', figure=_empty_fig(''),
                          config={'scrollZoom': True, 'displayModeBar': True,
                                  'modeBarButtonsToAdd': ['resetScale2d'],
                                  'modeBarButtonsToRemove': ['toImage', 'sendDataToCloud'],
                                  'displaylogo': False},
                          style={'height': '100%'}),
            ]),
        ]),
    ]),

    dcc.Store(id='selected-bin',    data=None),
    dcc.Store(id='selected-lonlat', data=None),
    dcc.Store(id='globe-scale',     data=1.0),
    dcc.Store(id='grid-res',        data=RES_FINE),
    dcc.Store(id='map-view',        data=VIEW_FLAT),
    dcc.Store(id='flat-proj',       data={'scale': 1.0, 'center_lat': 0.0, 'center_lon': 0.0}),
])

# ==============================================================================
# Callbacks
# ==============================================================================

@app.callback(
    Output('year-range-label', 'children'),
    Input('year-slider', 'value'),
)
def update_year_label(years):
    y0, y1 = years or [YEAR_MIN, YEAR_MAX]
    return f'{y0} – {y1}'


@app.callback(
    Output('globe',        'style'),
    Output('flatmap',      'style'),
    Output('zoom-box-wrap','style'),
    Output('flat-home-btn','style'),
    Output('map-view',     'data'),
    Input('view-radio', 'value'),
    prevent_initial_call=True,
)
def on_view_toggle(view):
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
    return globe_style, flat_style, zoom_style, home_style, view


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
    Output('globe',        'figure'),
    Output('flatmap',      'figure'),
    Output('grid-res',     'data'),
    Output('count-label',  'children'),
    Input('source-select', 'value'),
    Input('year-slider',   'value'),
    Input('res-radio',     'value'),
    State('selected-lonlat', 'data'),
    State('globe-scale',     'data'),
)
def on_map_controls(source_sel, years, res, sel_lonlat, scale):
    res = res or RES_FINE
    y0, y1 = years or [YEAR_MIN, YEAR_MAX]
    sel_lon = sel_lonlat['lon'] if sel_lonlat else None
    sel_lat = sel_lonlat['lat'] if sel_lonlat else None

    data  = query_profiles(source_sel, res, y0, y1)
    count = (f'showing {data["n_shown"]:,} of {data["n_total"]:,} profiles'
             if data else 'no profiles in this window')

    globe_fig = make_globe(data, source_sel, sel_lon, sel_lat, scale=scale or 1.0)
    flat_fig  = make_flat(data,  source_sel, sel_lon, sel_lat)
    return globe_fig, flat_fig, res, count


@app.callback(
    Output('selected-bin',    'data'),
    Output('selected-lonlat', 'data'),
    Output('globe',   'figure', allow_duplicate=True),
    Output('flatmap', 'figure', allow_duplicate=True),
    Input('globe',   'clickData'),
    Input('flatmap', 'clickData'),
    State('selected-lonlat', 'data'),
    prevent_initial_call=True,
)
def on_map_click(globe_click, flat_click, current_lonlat):
    from dash import ctx
    click_data = globe_click if ctx.triggered_id == 'globe' else flat_click
    if not click_data:
        return no_update, no_update, no_update, no_update
    points = click_data.get('points', [])
    if not points:
        return no_update, no_update, no_update, no_update
    p   = points[0]
    bid = p.get('customdata')
    if bid is None:
        return no_update, no_update, no_update, no_update

    sel_lon = p.get('lon')
    sel_lat = p.get('lat')
    sel_lonlat = {'lon': sel_lon, 'lat': sel_lat}

    sel_trace = _selection_trace(sel_lon, sel_lat)
    patched   = Patch()
    # Selection trace is always the last data trace (after source traces)
    # We rebuild it by updating via index -1 (last trace before coast on flat)
    # Simpler: just update both figs by patching last scatter trace.
    # Since trace count varies with sources present, use a full Patch on the last trace.
    patched['data'][-1]['lat'] = [sel_lat]
    patched['data'][-1]['lon'] = [sel_lon]
    patched['data'][-1]['marker'] = dict(size=14, symbol='star', color='white',
                                         line=dict(color='black', width=1))
    return bid, sel_lonlat, patched, patched


@app.callback(
    Output('zt-obs',    'figure'),
    Output('zt-anom',   'figure'),
    Output('bin-info',  'children'),
    Input('selected-bin',  'data'),
    Input('source-select', 'value'),
    Input('field-radio',   'value'),
    Input('anom-radio',    'value'),
    Input('grid-res',      'data'),
)
def update_plots(bin_id, source_sel, field, anom_ref, res):
    res = res or RES_FINE
    b_dict, c_dict = _grid(res)
    res_label = '10242-bin' if res == RES_FINE else '2562-bin'

    if bin_id is None or bin_id not in b_dict:
        empty = _empty_fig('Click a profile dot on the map')
        return empty, _empty_fig(''), ''

    resolved = resolve_bin(b_dict[bin_id], source_sel, res)
    lon, lat  = c_dict.get(bin_id, (np.nan, np.nan))
    n_months  = len(resolved['months'])

    if n_months == 0:
        src_note = 'all sources' if source_sel == SOURCE_ALL else source_sel
        msg = _empty_fig(f'Bin {bin_id} ({res_label}) has no {src_note} data')
        return msg, _empty_fig(''), html.Div(f'No {src_note} profiles in this bin.',
                                              style={'color': DIM, 'fontSize': '12px'})

    obs_fig  = make_zt_figure(resolved, field, depth)
    anom_var = f'{field}_minus_clim' if anom_ref == 'minus_clim' else f'{field}_anom'
    anom_fig = make_zt_figure(resolved, anom_var, depth)

    src_note = 'All sources' if source_sel == SOURCE_ALL else source_sel
    n_srcs   = len(b_dict[bin_id]['by_source'])
    info = [
        html.Div(f'Bin {bin_id} ({res_label})', style={'color': ACC, 'fontWeight': 'bold'}),
        html.Div(src_note, style={'color': DIM}),
        html.Div(f'lon {lon:.2f}°'),
        html.Div(f'lat {lat:.2f}°'),
        html.Div(f'{n_months} months with data'),
        html.Div(f'{resolved["n_profiles"]} profiles total'),
        html.Div(f'{n_srcs} source(s) in bin', style={'color': DIM}),
    ]
    return obs_fig, anom_fig, info


@app.callback(
    Output('zt-anom', 'figure', allow_duplicate=True),
    Input('zt-obs', 'relayoutData'),
    prevent_initial_call=True,
)
def sync_obs_to_anom(relayout):
    if not relayout:
        return no_update
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
    return patch if touched else no_update


@app.callback(
    Output('zt-obs', 'figure', allow_duplicate=True),
    Input('zt-anom', 'relayoutData'),
    prevent_initial_call=True,
)
def sync_anom_to_obs(relayout):
    if not relayout:
        return no_update
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
    return patch if touched else no_update


@app.callback(
    Output('flatmap',   'figure', allow_duplicate=True),
    Output('flat-proj', 'data'),
    Input('flatmap', 'relayoutData'),
    State('flat-proj', 'data'),
    prevent_initial_call=True,
)
def clamp_flat_zoom(relayout, proj_state):
    if not relayout:
        return no_update, no_update
    ps        = proj_state or {'scale': 1.0, 'center_lat': 0.0, 'center_lon': 0.0}
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
    Output('flatmap',   'figure', allow_duplicate=True),
    Output('flat-proj', 'data',   allow_duplicate=True),
    Input('flat-home-btn', 'n_clicks'),
    State('source-select', 'value'),
    State('year-slider',   'value'),
    State('grid-res',      'data'),
    State('selected-lonlat', 'data'),
    prevent_initial_call=True,
)
def flat_home(n_clicks, source_sel, years, res, sel_lonlat):
    y0, y1   = years or [YEAR_MIN, YEAR_MAX]
    res      = res or RES_FINE
    sel_lon  = sel_lonlat['lon'] if sel_lonlat else None
    sel_lat  = sel_lonlat['lat'] if sel_lonlat else None
    data     = query_profiles(source_sel or SOURCE_ALL, res, y0, y1)
    fig      = make_flat(data, source_sel or SOURCE_ALL, sel_lon, sel_lat,
                         uirevision=f'home-{n_clicks}')
    return fig, {'scale': 1.0, 'center_lat': 0.0, 'center_lon': 0.0}


if __name__ == '__main__':
    print("Starting app — open http://127.0.0.1:8050")
    app.run(debug=False, host='127.0.0.1', port=8050)
