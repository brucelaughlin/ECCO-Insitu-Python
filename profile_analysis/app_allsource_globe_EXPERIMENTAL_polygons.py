"""
Interactive multi-source globe app.

Generalizes app_pfl_globe.py to every NCEI-processed source (CTD_WOD, GLD_WOD,
MEOP, MRB_WOD, PFL, PFL_BGC, XBT_WOD, ...) via a Source selector. Reads the
partitioned pickle from build_allsource_bin_timeseries.py.

Run with:
    python app_allsource_globe.py
Then open http://127.0.0.1:8050 in a browser.

Controls:
  - Rotate the globe by dragging; zoom with the slider
  - Source selector: "All (combined)" or a single instrument
  - Click a bin to load its Z-T heatmaps and source table
  - Field T / S, anomaly reference (− climatology / − bin time-mean)
  - Probability threshold (>=1/2/3) and denominator basis

The pickle stores data PARTITIONED BY SOURCE. The combined ("All") view is
derived here on the fly by profile-weighted pooling of the per-source monthly
means (weight = #profiles that source contributed that month, applied per
depth level so a level missing in one source is filled by the others).

Requires: dash, plotly
Pickle produced by: build_allsource_bin_timeseries.py
"""

import pickle
import numpy as np
import plotly.graph_objects as go
from dash import Dash, dcc, html, Input, Output, State, dash_table, Patch, no_update
from pathlib import Path
from scipy.spatial import SphericalVoronoi

# ==============================================================================
# Config
# ==============================================================================

PICKLE_FILE = Path('/Users/brucel/ecco/yip/profile_data/z_profile_file_analysis/allsource_bin_timeseries.pkl')

SOURCE_ALL = 'ALL'   # sentinel for the combined view

# d3-geo (Plotly's geo engine) decides a spherical polygon's interior from ring
# winding, and its convention is the OPPOSITE of the GeoJSON RFC. If cells fill
# their complement (whole globe solid), flip this sign.
#   +1 : ring wound so area-vector points OUTWARD (GeoJSON/CCW-from-outside)
#   -1 : the opposite (what d3-geo usually wants)
CELL_WINDING = -1

# ==============================================================================
# Load data
# ==============================================================================

print("Loading pickle...", flush=True)
with open(PICKLE_FILE, 'rb') as f:
    data = pickle.load(f)

depth        = data['depth']           # (97,) float32
bins         = data['bins']            # dict: bin_id -> bin data (by_source)
bin_centres  = data['bin_centres']     # dict: bin_id -> (lon, lat)
SOURCES      = data['sources']         # sorted list of source names

print(f"  {len(bins):,} bins, {len(SOURCES)} sources: {', '.join(SOURCES)}")

DATASET_FIRST_YM = data.get('dataset_first_ym')
DATASET_LAST_YM  = data.get('dataset_last_ym')
DATASET_MONTHS   = data.get('total_months', 1)

N_DEPTH = len(depth)

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


print("Building geodesic cell polygons (SphericalVoronoi)...", flush=True)
CELL_GEOJSON, CELL_IDS = _build_cell_geojson(bin_centres)
print(f"  {len(CELL_IDS):,} cells built")


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


def _source_entries(month, source_sel):
    """List of by_source entries for the selected source (or all) in a month."""
    bs = month['by_source']
    if source_sel == SOURCE_ALL:
        return list(bs.values())
    entry = bs.get(source_sel)
    return [entry] if entry is not None else []


def resolve_bin(bin_data, source_sel):
    """Collapse a partitioned bin to a flat shape for the given source (or All).

    Returns a dict shaped like the old PFL pickle's bin so make_zt_figure and
    make_source_table work unchanged:
      {'bin_id', 'months': {ym: {'T','S','Tclim','Sclim','sources','n'}},
       'T_mean','S_mean','hit_counts','first_ym','last_ym', 'n_profiles'}
    """
    months = {}
    for ym, m in bin_data['months'].items():
        entries = _source_entries(m, source_sel)
        if not entries:
            continue
        # flatten source-profile records, tagging each with its source name
        recs = []
        for src_name, e in ((s, m['by_source'][s]) for s in m['by_source']
                            if source_sel in (SOURCE_ALL, s)):
            for r in e['sources']:
                rr = dict(r)
                rr['source'] = src_name
                recs.append(rr)
        months[ym] = {
            'T':      _pool_level_weighted(entries, 'T'),
            'S':      _pool_level_weighted(entries, 'S'),
            'Tclim':  _pool_level_weighted(entries, 'Tclim'),
            'Sclim':  _pool_level_weighted(entries, 'Sclim'),
            'sources': recs,
            'n':      sum(e['n'] for e in entries),
        }

    if source_sel == SOURCE_ALL:
        T_mean     = bin_data['T_mean']
        S_mean     = bin_data['S_mean']
        hit_counts = bin_data['hit_counts']
        first_ym   = bin_data['first_ym']
        last_ym    = bin_data['last_ym']
        n_profiles = sum(r['n_profiles'] for r in bin_data['by_source'].values())
    else:
        roll = bin_data['by_source'][source_sel]
        T_mean     = roll['T_mean']
        S_mean     = roll['S_mean']
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

def _visible_ids(source_sel):
    return [b for b in sorted(bins.keys())
            if b in bin_centres and bin_has_source(bins[b], source_sel)]


def _build_hover(ids_clean, source_sel, n_threshold, basis, prob_label):
    src_note = 'All sources' if source_sel == SOURCE_ALL else source_sel
    return [
        f'Bin {b} — {src_note}'
        f'<br>lon={bin_centres[b][0]:.2f}, lat={bin_centres[b][1]:.2f}'
        f'<br>P({prob_label} prof/month) = {bin_probability(bins[b], source_sel, n_threshold, basis):.3f}'
        for b in ids_clean
    ]


def _globe_arrays(source_sel, n_threshold, basis, selected_bin_id=None):
    """Per-visible-cell arrays for the choropleth: locations (bin ids), z
    (probability), and hover text."""
    ids_clean = _visible_ids(source_sel)
    probs     = [bin_probability(bins[b], source_sel, n_threshold, basis) for b in ids_clean]
    prob_label = next(k for k, v in PROB_KEYS.items() if v == n_threshold)
    hover     = _build_hover(ids_clean, source_sel, n_threshold, basis, prob_label)
    return dict(ids=ids_clean, probs=probs, hover=hover, prob_label=prob_label)


def _selected_trace(selected_bin_id):
    """A thin Scattergeo outline of the selected cell, drawn on top of the
    choropleth so the selection reads clearly. Empty when nothing is selected."""
    if selected_bin_id is None:
        return go.Scattergeo(lon=[], lat=[], mode='lines', hoverinfo='skip',
                             showlegend=False, name='selection')
    # find the cell ring from the prebuilt GeoJSON
    ring_lon, ring_lat = [], []
    for feat in CELL_GEOJSON['features']:
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
               n_threshold=1, basis=BASIS_DAY0, scale=1.0):
    a = _globe_arrays(source_sel, n_threshold, basis, selected_bin_id)

    fig = go.Figure(go.Choropleth(
        geojson=CELL_GEOJSON,
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
    fig.add_trace(_selected_trace(selected_bin_id))

    fig.update_geos(
        projection_type='orthographic',
        projection_scale=scale,
        projection_rotation=dict(lon=0, lat=20, roll=0),
        showland=True,   landcolor='#2a2a2a',
        showocean=True,  oceancolor='#0d1b2a',
        showcoastlines=True, coastlinecolor='#555',
        showframe=False,
        bgcolor='#1a1a1a',
    )
    fig.update_layout(
        margin=dict(l=0, r=0, t=0, b=0),
        paper_bgcolor='#1a1a1a',
        geo=dict(bgcolor='#1a1a1a'),
        autosize=True,
        dragmode=False,
        uirevision='globe',   # keeps rotation/zoom across figure updates
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
    lon, lat = bin_centres.get(bid, (np.nan, np.nan))
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

            // Wheel zoom: adjust ONLY geo.projection.scale, clamped, leaving
            // rotation/roll untouched. This deliberately bypasses Plotly's
            // built-in geo scrollZoom, which couples lon/lat/roll/scale into one
            // gesture and glitches over the poles (the original bug).
            const SCALE_MIN = 1.0, SCALE_MAX = 6.0;
            const WHEEL_STEP = 0.0015;   // per wheel delta unit; feels natural
            inner.addEventListener('wheel', function(e) {
                e.preventDefault();          // stop the page from scrolling
                e.stopPropagation();
                const geo = inner._fullLayout.geo;
                const cur = (geo && geo.projection && geo.projection.scale) || 1.0;
                // multiplicative so each notch is a constant zoom ratio
                let next = cur * Math.exp(-e.deltaY * WHEEL_STEP);
                next = Math.max(SCALE_MIN, Math.min(SCALE_MAX, next));
                Plotly.relayout(inner, {'geo.projection.scale': next});
                // Sync the slider handle without triggering on_zoom: that
                // callback listens to `drag_value`, so setting `value` here
                // moves the handle but does NOT loop back into a relayout.
                if (window.dash_clientside && window.dash_clientside.set_props) {
                    window.dash_clientside.set_props('zoom-slider',
                        {value: Math.round(next * 10) / 10});
                }
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

SCALE_MIN, SCALE_MAX = 1.0, 6.0

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
        html.Span('Pick a source, then click a bin to inspect its Z-T structure.',
                  style={'fontSize': '12px', 'color': '#888'}),
    ]),

    html.Div(style={'display': 'flex', 'gap': '8px', 'flex': '0 0 45vh',
                    'minHeight': '0'}, children=[

        html.Div(style={'flex': '1', 'backgroundColor': MID,
                        'borderRadius': '4px', 'padding': '4px',
                        'minHeight': '0', 'position': 'relative'}, children=[
            dcc.Graph(id='globe', figure=make_globe(),
                      config={'scrollZoom': False, 'doubleClick': False,
                              'displaylogo': False, 'displayModeBar': False},
                      style={'height': '100%'}),
            html.Div(className='zoom-box',
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
                html.Label('Source',
                           style={'fontSize': '12px', 'color': '#aaa',
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
                           style={'fontSize': '12px', 'color': '#aaa',
                                  'marginBottom': '6px', 'display': 'block'}),
                dcc.RadioItems(
                    id='prob-radio',
                    options=[{'label': f' {k}', 'value': k} for k in PROB_KEYS],
                    value='>=1',
                    labelStyle={'display': 'block', 'fontSize': '12px',
                                'marginBottom': '5px', 'cursor': 'pointer'},
                    inputStyle={'marginRight': '6px'},
                ),
            ]),

            html.Div([
                html.Label('Probability basis',
                           style={'fontSize': '12px', 'color': '#aaa',
                                  'marginBottom': '6px', 'display': 'block'}),
                dcc.RadioItems(
                    id='basis-radio',
                    options=[
                        {'label': ' Full dataset span', 'value': BASIS_DAY0},
                        {'label': ' Since bin\'s first obs', 'value': BASIS_SINCE_FIRST},
                    ],
                    value=BASIS_DAY0,
                    labelStyle={'display': 'block', 'fontSize': '12px',
                                'marginBottom': '5px', 'cursor': 'pointer'},
                    inputStyle={'marginRight': '6px'},
                ),
            ]),

            html.Div([
                html.Label('Field', style={'fontSize': '12px', 'color': '#aaa',
                                           'marginBottom': '6px', 'display': 'block'}),
                dcc.RadioItems(
                    id='field-radio',
                    options=[
                        {'label': ' Temperature (T)', 'value': 'T'},
                        {'label': ' Salinity (S)',    'value': 'S'},
                    ],
                    value='T',
                    labelStyle={'display': 'block', 'fontSize': '12px',
                                'marginBottom': '5px', 'cursor': 'pointer'},
                    inputStyle={'marginRight': '6px'},
                ),
            ]),

            html.Div([
                html.Label('Anomaly reference', style={'fontSize': '12px', 'color': '#aaa',
                                                       'marginBottom': '6px', 'display': 'block'}),
                dcc.RadioItems(
                    id='anom-radio',
                    options=[
                        {'label': ' − climatology',    'value': 'minus_clim'},
                        {'label': ' − bin time-mean',  'value': 'anom'},
                    ],
                    value='minus_clim',
                    labelStyle={'display': 'block', 'fontSize': '12px',
                                'marginBottom': '5px', 'cursor': 'pointer'},
                    inputStyle={'marginRight': '6px'},
                ),
            ]),

            html.Div([
                dcc.Checklist(
                    id='show-sources',
                    options=[{'label': ' Show source table', 'value': 'show'}],
                    value=[],
                    labelStyle={'fontSize': '12px', 'cursor': 'pointer'},
                    inputStyle={'marginRight': '6px'},
                ),
            ]),

            html.Div(id='bin-info', style={'fontSize': '11px', 'color': '#888',
                                            'lineHeight': '1.6'}),
        ]),
    ]),

    html.Div(style={'display': 'flex', 'gap': '8px', 'flex': '1',
                    'minHeight': '0', 'overflow': 'hidden'}, children=[

        html.Div(style={'flex': '1', 'backgroundColor': MID,
                        'borderRadius': '4px', 'padding': '4px',
                        'minWidth': '0'}, children=[
            dcc.Graph(id='zt-obs', figure=go.Figure(),
                      config={'scrollZoom': True, 'displayModeBar': True,
                              'modeBarButtonsToAdd': ['resetScale2d'],
                              'modeBarButtonsToRemove': ['toImage', 'sendDataToCloud'],
                              'displaylogo': False},
                      style={'height': '100%'}),
        ]),

        html.Div(id='zt-anom-container',
                 style={'flex': '1', 'backgroundColor': MID,
                        'borderRadius': '4px', 'padding': '4px',
                        'minWidth': '0'}, children=[
            dcc.Graph(id='zt-anom', figure=go.Figure(),
                      config={'scrollZoom': True, 'displayModeBar': True,
                              'modeBarButtonsToAdd': ['resetScale2d'],
                              'modeBarButtonsToRemove': ['toImage', 'sendDataToCloud'],
                              'displaylogo': False},
                      style={'height': '100%'}),
        ]),

        html.Div(id='source-container',
                 style={'display': 'none', 'flex': '1', 'backgroundColor': MID,
                        'borderRadius': '4px', 'padding': '8px',
                        'flexDirection': 'column', 'overflow': 'hidden',
                        'minWidth': '0'}, children=[
            html.Div('Source profiles', style={'fontSize': '12px', 'color': ACC,
                                                'marginBottom': '6px'}),
            html.Div(id='source-table-container', style={'flex': '1', 'overflowY': 'auto'}),
        ]),
    ]),

    dcc.Store(id='selected-bin', data=None),
    dcc.Store(id='globe-scale', data=1.0),
])

# ==============================================================================
# Callbacks
# ==============================================================================

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
    Output('globe', 'figure', allow_duplicate=True),
    Input('source-select', 'value'),
    Input('prob-radio', 'value'),
    Input('basis-radio', 'value'),
    State('selected-bin', 'data'),
    prevent_initial_call=True,
)
def on_globe_controls(source_sel, prob_label, basis, selected_bin):
    """Recolour/refilter the choropleth (which cells are visible + their z
    change with source), leaving layout (rotation/zoom) untouched."""
    n_threshold = PROB_KEYS[prob_label]
    a = _globe_arrays(source_sel, n_threshold, basis, selected_bin)
    patched = Patch()
    patched['data'][0]['locations']                           = a['ids']
    patched['data'][0]['z']                                   = a['probs']
    patched['data'][0]['customdata']                          = a['ids']
    patched['data'][0]['colorbar']['title']['text']           = f'P({prob_label})'
    patched['data'][0]['text']                                = a['hover']
    return patched


@app.callback(
    Output('selected-bin', 'data'),
    Output('globe', 'figure'),
    Input('globe', 'clickData'),
    State('selected-bin', 'data'),
    State('source-select', 'value'),
)
def on_globe_click(click_data, current_bin, source_sel):
    if click_data is None:
        return current_bin, no_update
    points = click_data.get('points', [])
    if not points:
        return current_bin, no_update
    # choropleth click returns the bin id via 'location' (and customdata)
    p = points[0]
    bid = p.get('customdata', p.get('location'))
    if bid is None or bid not in bins:
        return current_bin, no_update

    # Update only the selection-outline trace (data[1]); the choropleth
    # (data[0]) is untouched so no recolour fl: no flash.
    patched = Patch()
    sel = _selected_trace(bid)
    patched['data'][1]['lon'] = list(sel['lon'])
    patched['data'][1]['lat'] = list(sel['lat'])
    return bid, patched


def _empty_fig(message):
    fig = go.Figure()
    fig.update_layout(
        paper_bgcolor='#1a1a1a', plot_bgcolor='#111',
        xaxis=dict(visible=False), yaxis=dict(visible=False),
        annotations=[dict(text=message, xref='paper', yref='paper',
                          x=0.5, y=0.5, showarrow=False,
                          font=dict(size=14, color='#555'))],
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
)
def update_plots(bin_id, source_sel, field, anom_ref, prob_label, basis):
    if bin_id is None or bin_id not in bins:
        empty = _empty_fig('Click a bin on the globe')
        return (empty, _empty_fig(''),
                html.Div('No bin selected.', style={'color': '#555', 'fontSize': '12px'}),
                '')

    resolved = resolve_bin(bins[bin_id], source_sel)
    lon, lat = bin_centres.get(bin_id, (np.nan, np.nan))
    n_months = len(resolved['months'])
    n_profs  = resolved['n_profiles']

    if n_months == 0:
        src_note = 'all sources' if source_sel == SOURCE_ALL else source_sel
        msg = _empty_fig(f'Bin {bin_id} has no {src_note} data')
        return (msg, _empty_fig(''),
                html.Div(f'No {src_note} profiles in this bin.',
                         style={'color': '#555', 'fontSize': '12px'}),
                html.Div(f'Bin {bin_id}: no {src_note} data', style={'color': '#888'}))

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
    prob_val    = bin_probability(bins[bin_id], source_sel, n_threshold, basis)
    basis_note  = 'since first obs' if basis == BASIS_SINCE_FIRST else 'full span'
    src_note    = 'All sources' if source_sel == SOURCE_ALL else source_sel
    n_srcs      = len(bins[bin_id]['by_source'])
    info = [
        html.Div(f'Bin {bin_id}', style={'color': ACC, 'fontWeight': 'bold'}),
        html.Div(src_note, style={'color': '#aaa'}),
        html.Div(f'lon {lon:.2f}°'),
        html.Div(f'lat {lat:.2f}°'),
        html.Div(f'{n_months} months'),
        html.Div(f'{n_profs} profiles'),
        html.Div(f'{n_srcs} source(s) in bin', style={'color': '#666'}),
        html.Div(f'P({prob_glyph}) = {prob_val:.3f}'),
        html.Div(f'({basis_note})', style={'color': '#666'}),
    ]

    return obs_fig, anom_fig, table, info


if __name__ == '__main__':
    print("Starting app — open http://127.0.0.1:8050")
    app.run(debug=False, host='127.0.0.1', port=8050)
