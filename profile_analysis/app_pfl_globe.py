"""
Interactive PFL globe app.

Run with:
    python app_pfl_globe.py

Then open http://127.0.0.1:8050 in a browser.

Controls:
  - Rotate/zoom the globe by dragging
  - Click a bin to load its Z-T heatmap and source table
  - Radio buttons: T / S / T anomaly / S anomaly
  - Bin resolution toggle: 10242 bins (bin_key 'a')

Requires: dash, plotly  (pip install dash plotly)
Pickle produced by: build_pfl_bin_timeseries.py
"""

import pickle
import numpy as np
import pandas as pd
import plotly.graph_objects as go
import plotly.colors as pc
from dash import Dash, dcc, html, Input, Output, State, callback_context, dash_table, Patch, no_update
import json
from pathlib import Path

# ==============================================================================
# Config
# ==============================================================================

PICKLE_FILE = Path('/Users/brucel/ecco/yip/profile_data/z_profile_file_analysis/pfl_bin_timeseries.pkl')

# ==============================================================================
# Load data
# ==============================================================================

print("Loading pickle...", flush=True)
with open(PICKLE_FILE, 'rb') as f:
    data = pickle.load(f)

depth        = data['depth']           # (97,)
bins         = data['bins']            # dict: bin_id -> bin data
bin_centres  = data['bin_centres']     # dict: bin_id -> (lon, lat)

print(f"  {len(bins):,} bins loaded")

DATASET_FIRST_YM = data.get('dataset_first_ym')
DATASET_LAST_YM  = data.get('dataset_last_ym')
DATASET_MONTHS   = data.get('total_months', 1)


def _months_between(ym0, ym1):
    """Inclusive count of calendar months from ym0 to ym1, each (year, month)."""
    return (ym1[0] - ym0[0]) * 12 + (ym1[1] - ym0[1]) + 1


# Backfill ingredients if loading an older pickle without hit_counts/first_ym
for b in bins.values():
    if 'hit_counts' not in b:
        hc = {1: 0, 2: 0, 3: 0}
        for entry in b['months'].values():
            n = len(entry['sources'])
            if n >= 1: hc[1] += 1
            if n >= 2: hc[2] += 1
            if n >= 3: hc[3] += 1
        b['hit_counts'] = hc
        yms = sorted(b['months'].keys())
        b['first_ym'] = yms[0] if yms else None
        b['last_ym']  = yms[-1] if yms else None

# Threshold selector maps label -> N
PROB_KEYS = {'>=1': 1, '>=2': 2, '>=3': 3}

# Denominator-basis selector
BASIS_DAY0        = 'day0'          # full dataset span (same denom for all bins)
BASIS_SINCE_FIRST = 'since_first'   # bin's first obs -> dataset end


def bin_probability(b, n_threshold, basis):
    """Probability that this bin received >= n_threshold profiles in a month,
    under the chosen denominator basis."""
    hits = b['hit_counts'].get(n_threshold, 0)
    if basis == BASIS_SINCE_FIRST and b.get('first_ym') and DATASET_LAST_YM:
        denom = _months_between(b['first_ym'], DATASET_LAST_YM)
    else:
        denom = DATASET_MONTHS
    return hits / denom if denom > 0 else 0.0

# ==============================================================================
# Build globe figure
# ==============================================================================

def _build_hover(ids_clean, n_threshold, basis, prob_label):
    return [
        f'Bin {b}<br>lon={bin_centres[b][0]:.2f}, lat={bin_centres[b][1]:.2f}'
        f'<br>P({prob_label} prof/month) = {bin_probability(bins[b], n_threshold, basis):.3f}'
        f'<br>{len(bins[b]["months"])} months with data'
        for b in ids_clean
    ]


def make_globe(selected_bin_id=None, n_threshold=1, basis=BASIS_DAY0, scale=1.0):
    bin_ids   = sorted(bins.keys())
    lons      = [bin_centres[b][0] for b in bin_ids if b in bin_centres]
    lats      = [bin_centres[b][1] for b in bin_ids if b in bin_centres]
    ids_clean = [b for b in bin_ids if b in bin_centres]
    probs     = [bin_probability(bins[b], n_threshold, basis) for b in ids_clean]

    # Marker sizes proportional to probability, min visible size
    sizes = [max(3, p * 14) for p in probs]

    # Color by probability
    colorscale = 'RdYlGn'
    marker_colors = probs

    # Highlight selected bin
    outline_colors = ['rgba(255,255,255,0.9)' if b == selected_bin_id
                      else 'rgba(0,0,0,0)' for b in ids_clean]
    outline_widths = [2 if b == selected_bin_id else 0 for b in ids_clean]

    prob_label = next(k for k, v in PROB_KEYS.items() if v == n_threshold)
    hover = _build_hover(ids_clean, n_threshold, basis, prob_label)

    fig = go.Figure(go.Scattergeo(
        lon=lons,
        lat=lats,
        mode='markers',
        marker=dict(
            size=sizes,
            color=marker_colors,
            colorscale=colorscale,
            cmin=0, cmax=1,
            colorbar=dict(title=f'P({prob_label})', thickness=14, len=0.6),
            line=dict(color=outline_colors, width=outline_widths),
        ),
        text=hover,
        hoverinfo='text',
        customdata=ids_clean,
    ))

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
        dragmode=False,      # disable Plotly's built-in geo drag; we drive rotation ourselves
        uirevision='globe',  # keeps rotation/zoom across figure updates
    )
    return fig


# ==============================================================================
# Z-T heatmap helpers
# ==============================================================================

def make_zt_figure(bin_data, var, depth):
    """
    var: 'T' | 'S' | 'Tclim' | 'Sclim' | 'T_minus_clim' | 'S_minus_clim'
         | 'T_anom' | 'S_anom'
    Returns a plotly Figure.
    """
    months_dict = bin_data['months']
    sorted_yms  = sorted(months_dict.keys())
    if not sorted_yms:
        return go.Figure()

    x = [f'{ym[0]}-{ym[1]:02d}-15' for ym in sorted_yms]

    is_diff_clim = var in ('T_minus_clim', 'S_minus_clim')
    is_anom      = var in ('T_anom', 'S_anom')
    base_letter  = var[0]   # 'T' or 'S'

    if is_diff_clim:
        obs_key  = base_letter          # 'T' or 'S'
        clim_key = base_letter + 'clim' # 'Tclim' or 'Sclim'
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
        # Raw: 'T', 'S', 'Tclim', 'Sclim'
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
        x=x,
        y=depth,
        z=matrix,
        colorscale=colorscale,
        zmin=zmin, zmax=zmax,
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
                'Type':       s['type'],
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
app.title = 'PFL Argo Globe'

# Brighten the dcc.Slider track/handle so the zoom bar is visible on the dark globe
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

# Custom globe rotation: we drive the orthographic projection ourselves instead
# of using Plotly's built-in geo drag (which couples lon/lat/roll/zoom into one
# gesture and glitches over the poles). This handler:
#   - rotates longitude freely (wraps naturally)
#   - rotates latitude, hard-clamped to +/-LAT_CAP during the drag
#   - never touches roll (stays 0 -> globe can't flip/blank)
#   - never touches scale (zoom only changes via the slider)
#   - still lets genuine clicks through (small movement) for bin selection
app.clientside_callback(
    """
    function(_) {
        const LAT_CAP    = 85;    // deg; last few degrees near pole are unstable
        const DEG_PER_PX = 0.4;   // rotation sensitivity
        const CLICK_TOL  = 5;     // px; movement below this counts as a click

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

                // Horizontal drag -> longitude; vertical drag -> latitude.
                let newLon = startLon - dx * DEG_PER_PX;
                let newLat = startLat + dy * DEG_PER_PX;
                // Wrap longitude to [-180,180]; clamp latitude.
                newLon = ((newLon + 180) % 360 + 360) % 360 - 180;
                newLat = Math.max(-LAT_CAP, Math.min(LAT_CAP, newLat));

                Plotly.relayout(inner, {
                    'geo.projection.rotation.lon': newLon,
                    'geo.projection.rotation.lat': newLat,
                    'geo.projection.rotation.roll': 0,
                });
            }, true);

            window.addEventListener('mouseup', function() { dragging = false; }, true);

            // Suppress Plotly's click (bin select) if this was actually a drag.
            inner.addEventListener('click', function(e) {
                if (moved) { e.stopImmediatePropagation(); }
            }, true);

            // Belt-and-suspenders: kill double-click zoom.
            inner.on('plotly_doubleclick', function() { return false; });
        }
        setup();
        return window.dash_clientside.no_update;
    }
    """,
    Output('globe', 'id'),   # dummy output — we just need a trigger
    Input('globe', 'id'),
)

DARK = '#1a1a1a'
MID  = '#252525'
TEXT = '#dddddd'
ACC  = '#88ccff'

# Globe zoom scale bounds
SCALE_MIN, SCALE_MAX = 1.0, 6.0

app.layout = html.Div(style={'backgroundColor': DARK, 'color': TEXT,
                              'fontFamily': 'sans-serif', 'height': '100vh',
                              'display': 'flex', 'flexDirection': 'column',
                              'padding': '8px', 'gap': '8px'}, children=[

    # Title bar
    html.Div(style={'backgroundColor': MID, 'padding': '8px 14px',
                    'borderRadius': '4px', 'display': 'flex',
                    'alignItems': 'center', 'gap': '20px'}, children=[
        html.H2('Argo PFL — Interactive Bin Explorer',
                style={'margin': 0, 'fontSize': '16px', 'color': ACC}),
        html.Span('Click a bin on the globe to inspect its Z-T structure.',
                  style={'fontSize': '12px', 'color': '#888'}),
    ]),

    # Globe + controls row
    html.Div(style={'display': 'flex', 'gap': '8px', 'flex': '0 0 45vh',
                    'minHeight': '0'}, children=[

        # Globe
        html.Div(style={'flex': '1', 'backgroundColor': MID,
                        'borderRadius': '4px', 'padding': '4px',
                        'minHeight': '0', 'position': 'relative'}, children=[
            dcc.Graph(id='globe', figure=make_globe(),
                      config={'scrollZoom': False, 'doubleClick': False,
                              'displaylogo': False, 'displayModeBar': False},
                      style={'height': '100%'}),
            # Zoom slider — overlaid top-left of the globe, on a visible backdrop
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

        # Right-side controls
        html.Div(style={'width': '200px', 'backgroundColor': MID,
                        'borderRadius': '4px', 'padding': '14px',
                        'display': 'flex', 'flexDirection': 'column',
                        'gap': '14px', 'overflowY': 'auto', 'minHeight': '0'},
                 children=[

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

    # Lower row: observed Z-T (left) + anomaly Z-T OR source table (right)
    html.Div(style={'display': 'flex', 'gap': '8px', 'flex': '1',
                    'minHeight': '0', 'overflow': 'hidden'}, children=[

        # Left: observed Z-T
        html.Div(style={'flex': '1', 'backgroundColor': MID,
                        'borderRadius': '4px', 'padding': '4px',
                        'minWidth': '0'}, children=[
            dcc.Graph(id='zt-obs', figure=go.Figure(),
                      config={'scrollZoom': True,
                              'displayModeBar': True,
                              'modeBarButtonsToAdd': ['resetScale2d'],
                              'modeBarButtonsToRemove': ['toImage', 'sendDataToCloud'],
                              'displaylogo': False},
                      style={'height': '100%'}),
        ]),

        # Right: anomaly Z-T (default) — hidden when source table is shown
        html.Div(id='zt-anom-container',
                 style={'flex': '1', 'backgroundColor': MID,
                        'borderRadius': '4px', 'padding': '4px',
                        'minWidth': '0'}, children=[
            dcc.Graph(id='zt-anom', figure=go.Figure(),
                      config={'scrollZoom': True,
                              'displayModeBar': True,
                              'modeBarButtonsToAdd': ['resetScale2d'],
                              'modeBarButtonsToRemove': ['toImage', 'sendDataToCloud'],
                              'displaylogo': False},
                      style={'height': '100%'}),
        ]),

        # Right (alternate): source table — hidden by default
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

    # Store the currently selected bin_id and the globe zoom scale
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
    # Patch only the projection scale — rotation (geo.projection.rotation) untouched
    patched = Patch()
    patched['layout']['geo']['projection']['scale'] = scale
    return patched, scale


@app.callback(
    Output('globe', 'figure', allow_duplicate=True),
    Input('prob-radio', 'value'),
    Input('basis-radio', 'value'),
    prevent_initial_call=True,
)
def on_prob_change(prob_label, basis):
    n_threshold = PROB_KEYS[prob_label]
    ids_clean = [b for b in sorted(bins.keys()) if b in bin_centres]
    probs     = [bin_probability(bins[b], n_threshold, basis) for b in ids_clean]
    sizes     = [max(3, p * 14) for p in probs]

    patched = Patch()
    patched['data'][0]['marker']['color']          = probs
    patched['data'][0]['marker']['size']           = sizes
    patched['data'][0]['marker']['colorbar']['title']['text'] = f'P({prob_label})'
    patched['data'][0]['text']                     = _build_hover(ids_clean, n_threshold, basis, prob_label)
    return patched


@app.callback(
    Output('selected-bin', 'data'),
    Output('globe', 'figure'),
    Input('globe', 'clickData'),
    State('selected-bin', 'data'),
)
def on_globe_click(click_data, current_bin):
    if click_data is None:
        return current_bin, no_update

    points = click_data.get('points', [])
    if not points:
        return current_bin, no_update

    bid = points[0].get('customdata')
    if bid is None or bid not in bins:
        return current_bin, no_update

    # Use Patch to only update marker outlines — never touch layout/geo/zoom
    ids_clean = [b for b in sorted(bins.keys()) if b in bin_centres]
    patched = Patch()
    patched['data'][0]['marker']['line']['color'] = [
        'rgba(255,255,255,0.9)' if b == bid else 'rgba(0,0,0,0)'
        for b in ids_clean
    ]
    patched['data'][0]['marker']['line']['width'] = [
        2 if b == bid else 0
        for b in ids_clean
    ]
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
    """Translate a relayoutData dict into a Patch that mirrors x/y axis ranges."""
    if not relayout:
        return None
    patch = Patch()
    touched = False

    # Explicit range drag/zoom
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

    # Reset / autorange (double-click or reset-axes button)
    if relayout.get('xaxis.autorange'):
        patch['layout']['xaxis']['autorange'] = True
        touched = True
    if relayout.get('yaxis.autorange'):
        # y is reversed depth axis
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
    Input('field-radio', 'value'),
    Input('anom-radio', 'value'),
    Input('prob-radio', 'value'),
    Input('basis-radio', 'value'),
)
def update_plots(bin_id, field, anom_ref, prob_label, basis):
    if bin_id is None or bin_id not in bins:
        empty = _empty_fig('Click a bin on the globe')
        return (empty, _empty_fig(''),
                html.Div('No bin selected.', style={'color': '#555', 'fontSize': '12px'}),
                '')

    bin_data = bins[bin_id]
    lon, lat = bin_centres.get(bin_id, (np.nan, np.nan))
    n_months = len(bin_data['months'])
    n_profs  = sum(len(v['sources']) for v in bin_data['months'].values())

    # Left: observed field
    obs_fig = make_zt_figure(bin_data, field, depth)

    # Right: anomaly (relative to climatology or bin time-mean)
    anom_var = f'{field}_minus_clim' if anom_ref == 'minus_clim' else f'{field}_anom'
    anom_fig = make_zt_figure(bin_data, anom_var, depth)

    # Source table
    rows = make_source_table(bin_data)
    table = dash_table.DataTable(
        data=rows,
        columns=[{'name': c, 'id': c} for c in
                 ('Year-Month', 'Date', 'Type', 'Lon', 'Lat', 'File', 'Prof index')],
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
    prob_val    = bin_probability(bin_data, n_threshold, basis)
    basis_note  = 'since first obs' if basis == BASIS_SINCE_FIRST else 'full span'
    info = [
        html.Div(f'Bin {bin_id}', style={'color': ACC, 'fontWeight': 'bold'}),
        html.Div(f'lon {lon:.2f}°'),
        html.Div(f'lat {lat:.2f}°'),
        html.Div(f'{n_months} months'),
        html.Div(f'{n_profs} profiles'),
        html.Div(f'P({prob_glyph}) = {prob_val:.3f}'),
        html.Div(f'({basis_note})', style={'color': '#666'}),
    ]

    return obs_fig, anom_fig, table, info


if __name__ == '__main__':
    print("Starting app — open http://127.0.0.1:8050")
    app.run(debug=False, host='127.0.0.1', port=8050)
