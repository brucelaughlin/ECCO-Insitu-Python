"""
Map pieces shared by the profile map app (app_profile_map.py): the three flat-map
projections, the geodesic bin cells, the coastline, and the zoom/pan bookkeeping.

Projections
    world     natural earth, the whole globe
    arctic    polar stereographic centered on the North Pole, cut at POLAR_LAT_LIMIT N
    southern  polar stereographic centered on the South Pole, cut at POLAR_LAT_LIMIT S
"""

import numpy as np
import plotly.graph_objects as go
import cartopy.feature as cfeature
from scipy.spatial import SphericalVoronoi

PROJ_WORLD, PROJ_ARCTIC, PROJ_SOUTHERN = 'world', 'arctic', 'southern'
POLAR_LAT_LIMIT = 50          # polar maps open showing poleward of this latitude

# Plotly's stereographic map at scale 1 reaches the equator (radius ∝ tan(colat/2)),
# and lataxis_range does not crop it; so the polar maps open zoomed in until the
# POLAR_LAT_LIMIT circle spans ~90% of the map height (scale ≈ 2.5 for 50°).
POLAR_BASE_SCALE = 0.9 / np.tan(np.radians(90 - POLAR_LAT_LIMIT) / 2)

MAX_ZOOM = 40                 # most zoom-in, relative to the opening view

# minscale / maxscale make the browser itself stop the zoom at the opening view and at
# MAX_ZOOM (correcting it from the server fought the scroll wheel: the map bounced).
PROJECTIONS = {
    PROJ_WORLD:    dict(label='World', base_scale=1.0,
                        geo=dict(projection_type='natural earth', projection_scale=1.0,
                                 projection_minscale=1.0, projection_maxscale=MAX_ZOOM,
                                 lataxis_range=[-90, 90], lonaxis_range=[-180, 180])),
    PROJ_ARCTIC:   dict(label='Arctic', base_scale=POLAR_BASE_SCALE,
                        geo=dict(projection_type='stereographic', projection_scale=POLAR_BASE_SCALE,
                                 projection_minscale=POLAR_BASE_SCALE,
                                 projection_maxscale=POLAR_BASE_SCALE * MAX_ZOOM,
                                 projection_rotation=dict(lon=0, lat=90, roll=0))),
    PROJ_SOUTHERN: dict(label='Southern Ocean', base_scale=POLAR_BASE_SCALE,
                        geo=dict(projection_type='stereographic', projection_scale=POLAR_BASE_SCALE,
                                 projection_minscale=POLAR_BASE_SCALE,
                                 projection_maxscale=POLAR_BASE_SCALE * MAX_ZOOM,
                                 projection_rotation=dict(lon=0, lat=-90, roll=0))),
}

GEO_STYLE = dict(
    showland=True,    landcolor='#2a2a2a',
    showocean=True,   oceancolor='#0d1b2a',
    showlakes=True,   lakecolor='#3a3a3a',
    showrivers=True,  rivercolor='#3a3a3a',
    showframe=False,  bgcolor='#1a1a1a',
    showcoastlines=False,          # drawn by coastline_trace() instead (finer)
)

BG = '#1a1a1a'


def apply_projection(fig, proj, uirevision):
    """Set the geo layout of fig for projection key `proj`."""
    fig.update_geos(**PROJECTIONS[proj]['geo'], **GEO_STYLE)
    fig.update_layout(margin=dict(l=0, r=0, t=0, b=0), paper_bgcolor=BG,
                      autosize=True, dragmode='pan', uirevision=uirevision)
    return fig


# ------------------------------------------------------------------------------
# Geodesic bin cells
# ------------------------------------------------------------------------------

CELL_WINDING = -1


def build_cells(centers):
    """Spherical Voronoi cells of the bin centers.

    Returns (geojson FeatureCollection, ids, rings) where rings[bid] is the closed
    ring as an (n, 2) lon/lat array, longitudes unwrapped around the bin center.
    """
    ids = sorted(centers.keys())
    lon = np.radians([centers[b][0] for b in ids])
    lat = np.radians([centers[b][1] for b in ids])
    xyz = np.stack([np.cos(lat) * np.cos(lon), np.cos(lat) * np.sin(lon), np.sin(lat)], axis=-1)
    sv = SphericalVoronoi(xyz, radius=1.0, center=np.zeros(3))
    sv.sort_vertices_of_regions()

    features, rings = [], {}
    for k, bid in enumerate(ids):
        region = sv.regions[k]
        verts3d = sv.vertices[region]
        area_vec = np.zeros(3)
        for i in range(len(verts3d)):
            area_vec += np.cross(verts3d[i], verts3d[(i + 1) % len(verts3d)])
        if np.dot(area_vec, xyz[k]) * CELL_WINDING < 0:
            region = region[::-1]
        v = sv.vertices[region]
        lons = np.degrees(np.arctan2(v[:, 1], v[:, 0]))
        lats = np.degrees(np.arcsin(np.clip(v[:, 2], -1.0, 1.0)))
        c_lon = centers[bid][0]
        lons = c_lon + ((lons - c_lon + 180.0) % 360.0 - 180.0)
        ring = np.column_stack([lons, lats])
        ring = np.vstack([ring, ring[:1]])
        rings[bid] = ring
        features.append({'type': 'Feature', 'id': int(bid),            # 0.001° ≈ 100 m: halves the figure
                         'geometry': {'type': 'Polygon', 'coordinates': [np.round(ring, 3).tolist()]},
                         'properties': {}})
    return {'type': 'FeatureCollection', 'features': features}, ids, rings


def outline_trace(rings, keep, color='rgba(200,200,200,0.22)', width=0.5, name='_bins'):
    """Faint outlines of the bins in `keep` as one lines trace (rings broken at
    antimeridian jumps)."""
    out_lon, out_lat = [], []
    for bid in sorted(keep):
        ring = rings.get(bid)
        if ring is None:
            continue
        r_lon = (ring[:, 0] + 180.0) % 360.0 - 180.0
        for i in range(len(ring)):
            if i and abs(r_lon[i] - r_lon[i - 1]) > 180:
                out_lon.append(None); out_lat.append(None)
            out_lon.append(float(r_lon[i])); out_lat.append(float(ring[i, 1]))
        out_lon.append(None); out_lat.append(None)
    return go.Scattergeo(lon=out_lon, lat=out_lat, mode='lines',
                         line=dict(color=color, width=width),
                         hoverinfo='skip', showlegend=False, name=name)


def coastline_trace():
    lons, lats = [], []
    for g in cfeature.NaturalEarthFeature('physical', 'coastline', '50m').geometries():
        parts = [g] if g.geom_type == 'LineString' else list(g.geoms)
        for part in parts:
            x, y = part.xy
            lons.extend(list(x) + [None]); lats.extend(list(y) + [None])
    return go.Scattergeo(lon=lons, lat=lats, mode='lines',
                         line=dict(color='#4a9eff', width=0.9),
                         hoverinfo='skip', showlegend=False, name='_coast')


# ------------------------------------------------------------------------------
# Zoom / pan bookkeeping
# ------------------------------------------------------------------------------

def default_view(proj):
    """View state: zoom scale and the lon/lat at the center of the map."""
    lat = {PROJ_ARCTIC: 90.0, PROJ_SOUTHERN: -90.0}.get(proj, 0.0)
    return {'proj': proj, 'scale': PROJECTIONS[proj]['base_scale'], 'lon': 0.0, 'lat': lat}


def zoom_factor(state):
    """How far the view is zoomed in relative to the projection's opening view."""
    return state['scale'] / PROJECTIONS[state['proj']]['base_scale']


def update_view(relayout, state, proj):
    """Fold a map relayoutData event into the view state.

    Plotly reports pan as geo.center.* on the world map and as
    geo.projection.rotation.* on the (clipped) polar maps; both are read.
    Returns (new_state, clamp) where clamp holds values to patch back into the
    layout ('center_lat'): the world map panned past the poles. Else {}. Polar maps
    pan freely. (Zoom limits are minscale / maxscale, in the browser.)
    """
    st = dict(state) if state and state.get('proj') == proj else default_view(proj)
    if not relayout:
        return st, {}
    if relayout.get('geo.projection.scale') is not None:
        st['scale'] = float(relayout['geo.projection.scale'])
    for k in ('geo.center.lon', 'geo.projection.rotation.lon'):
        if relayout.get(k) is not None:
            st['lon'] = float(relayout[k])
    for k in ('geo.center.lat', 'geo.projection.rotation.lat'):
        if relayout.get(k) is not None:
            st['lat'] = float(relayout[k])
    clamp = {}
    base = PROJECTIONS[proj]['base_scale']
    st['scale'] = max(st['scale'], base)
    if proj == PROJ_WORLD:
        max_lat = max(0.0, 81.0 - 81.0 * base / st['scale'])
        if abs(st['lat']) > max_lat:
            st['lat'] = float(np.clip(st['lat'], -max_lat, max_lat))
            clamp['center_lat'] = st['lat']
    return st, clamp


def zoom_bucket(state):
    """Zoom rounded down to 1, 2, 4 or 8: things that scale with zoom (edge widths)
    change in steps, so the map is not redrawn on every scroll tick."""
    return int(2 ** np.clip(np.floor(np.log2(max(zoom_factor(state), 1.0))), 0, 3)) if state else 1
