"""
Sub-polygons for the profile map app: split one geodesic bin into n cells, and draw
a bin's outline slightly inside it.

The subdivision is the one from the original matplotlib tool
(geodesic_binning/binning/geodesic_binning_utilities_binning.py, _fill_polygon_voronoi
and helpers): one seed per horizontal band (so cells stack south to north), two
Lloyd relaxation passes, then Voronoi cells clipped to the bin. Changes from the
original:
  * the work is done in a gnomonic (tangent-plane) projection centered on the bin,
    so bins at the poles and across the antimeridian come out right (the original
    worked in plain lon/lat and skipped or distorted them); bin edges and cell
    edges are great circles, which are straight lines in this projection;
  * the sampling loops are bounded (the original could loop forever on a
    degenerate polygon);
  * a cell that fails is dropped, rather than replaced by the whole bin.

All rings are closed (n, 2) arrays of lon / lat in degrees, longitudes unwrapped
around the bin center, wound the same way as the bin ring they came from.
"""

import numpy as np
import shapely
from matplotlib.path import Path as MplPath
from scipy.spatial import Voronoi

_MAX_BATCHES = 50          # rejection-sampling attempts before falling back


# ------------------------------------------------------------------------------
# Gnomonic projection around (lon0, lat0)
# ------------------------------------------------------------------------------

def _to_plane(lonlat, lon0, lat0):
    lam, phi = np.radians(lonlat[:, 0]), np.radians(lonlat[:, 1])
    lam0, phi0 = np.radians(lon0), np.radians(lat0)
    dlam = lam - lam0
    cos_c = np.sin(phi0) * np.sin(phi) + np.cos(phi0) * np.cos(phi) * np.cos(dlam)
    x = np.cos(phi) * np.sin(dlam) / cos_c
    y = (np.cos(phi0) * np.sin(phi) - np.sin(phi0) * np.cos(phi) * np.cos(dlam)) / cos_c
    return np.column_stack([x, y])


def _to_lonlat(xy, lon0, lat0):
    x, y = xy[:, 0], xy[:, 1]
    lam0, phi0 = np.radians(lon0), np.radians(lat0)
    rho = np.hypot(x, y)
    c = np.arctan(rho)
    with np.errstate(invalid='ignore', divide='ignore'):
        phi = np.where(rho > 0,
                       np.arcsin(np.cos(c) * np.sin(phi0) + y * np.sin(c) * np.cos(phi0) / rho),
                       phi0)
    lam = lam0 + np.arctan2(x * np.sin(c), rho * np.cos(phi0) * np.cos(c) - y * np.sin(phi0) * np.sin(c))
    lon = np.degrees(lam)
    lon = lon0 + ((lon - lon0 + 180.0) % 360.0 - 180.0)          # unwrap around the center
    return np.column_stack([lon, np.degrees(phi)])


def _signed_area(xy):
    x, y = xy[:, 0], xy[:, 1]
    return 0.5 * float(np.sum(x[:-1] * y[1:] - x[1:] * y[:-1]))


def _closed(xy):
    return xy if np.allclose(xy[0], xy[-1]) else np.vstack([xy, xy[:1]])


# ------------------------------------------------------------------------------
# Planar subdivision (lifted from the original tool, loops bounded)
# ------------------------------------------------------------------------------

def _sample_points(path, bounds, n, rng):
    """Rejection sampling of n interior points (fewer if the polygon is degenerate)."""
    minx, miny, maxx, maxy = bounds
    pts = []
    for _ in range(_MAX_BATCHES):
        batch = rng.uniform([minx, miny], [maxx, maxy], size=(max(4 * n, 64), 2))
        pts.extend(batch[path.contains_points(batch)].tolist())
        if len(pts) >= n:
            break
    return np.array(pts[:n]).reshape(-1, 2)


def _seeds_y_stratified(path, bounds, n, rng, fallback):
    """One seed per horizontal band, so cells stack south to north (key for small n)."""
    minx, miny, maxx, maxy = bounds
    band = (maxy - miny) / n
    seeds = []
    for i in range(n):
        lo = miny + i * band
        pt = None
        for _ in range(_MAX_BATCHES):
            batch = rng.uniform([minx, lo], [maxx, lo + band], size=(64, 2))
            inside = path.contains_points(batch)
            if inside.any():
                pt = batch[inside][0]
                break
        if pt is None:                       # band outside a concave polygon: anywhere inside
            extra = _sample_points(path, bounds, 1, rng)
            pt = extra[0] if len(extra) else fallback
        seeds.append(pt)
    return np.array(seeds)


def _voronoi_cells(seeds, poly):
    """Voronoi cells of seeds clipped to poly (shapely); None where a cell fails."""
    minx, miny, maxx, maxy = poly.bounds
    m = max(maxx - minx, maxy - miny)
    ghosts = np.array([[minx - m, miny - m], [maxx + m, miny - m],
                       [minx - m, maxy + m], [maxx + m, maxy + m]])
    vor = Voronoi(np.vstack([seeds, ghosts]))
    regions = []
    for i in range(len(seeds)):
        reg = vor.regions[vor.point_region[i]]
        regions.append(None if (-1 in reg or len(reg) < 3) else shapely.Polygon(vor.vertices[reg]))
    ok = [i for i, r in enumerate(regions) if r is not None]
    clipped = shapely.intersection(np.array([regions[i] for i in ok], dtype=object), poly) if ok else []
    cells = [None] * len(seeds)
    for i, g in zip(ok, clipped):
        if g.is_empty:
            continue
        if g.geom_type == 'MultiPolygon' or g.geom_type == 'GeometryCollection':
            parts = [p for p in getattr(g, 'geoms', []) if p.geom_type == 'Polygon']
            if not parts:
                continue
            g = max(parts, key=lambda p: p.area)        # keep the largest piece
        if g.geom_type == 'Polygon' and g.area > 0:
            cells[i] = g
    return cells


def _lloyd_step(seeds, path, poly, rng):
    """One Lloyd pass: move each seed to its cell's centroid (inside the polygon)."""
    cells = _voronoi_cells(seeds, poly)
    new = seeds.copy()
    for i, cell in enumerate(cells):
        if cell is None:
            continue
        c = np.array([cell.centroid.x, cell.centroid.y])
        if path.contains_points([c])[0]:
            new[i] = c
        else:                                             # concave polygon: a random interior point
            extra = _sample_points(path, poly.bounds, 1, rng)
            if len(extra):
                new[i] = extra[0]
    return new


def _subdivide_plane(xy, n, seed, lloyd_iters=2):
    poly = shapely.Polygon(xy)
    if not poly.is_valid or poly.area <= 0:
        return []
    if n <= 1:
        return [xy]
    path = MplPath(xy)
    rng = np.random.default_rng(seed)
    rep = np.array(poly.representative_point().coords[0])
    seeds = _seeds_y_stratified(path, poly.bounds, n, rng, rep)
    for _ in range(lloyd_iters):
        seeds = _lloyd_step(seeds, path, poly, rng)
    return [np.asarray(c.exterior.coords) for c in _voronoi_cells(seeds, poly) if c is not None]


# ------------------------------------------------------------------------------
# Public
# ------------------------------------------------------------------------------

def subdivide(ring, center, n, seed):
    """Split a bin into (up to) n cells.

    ring   : closed (k, 2) lon/lat ring of the bin
    center : (lon, lat) of the bin center (the projection center)
    n      : number of cells
    seed   : random seed (use the bin id, so a bin always looks the same)
    Returns a list of closed lon/lat rings (fewer than n if a cell failed), wound
    like `ring`.
    """
    lon0, lat0 = center
    xy = _closed(_to_plane(np.asarray(ring, dtype=float), lon0, lat0))
    sign = np.sign(_signed_area(xy))
    out = []
    for cell in _subdivide_plane(xy, int(n), int(seed)):
        cell = _closed(cell)
        if np.sign(_signed_area(cell)) != sign:
            cell = cell[::-1]
        out.append(_to_lonlat(cell, lon0, lat0))
    return out


def order_south_to_north(rings):
    """Indices of rings sorted by centroid latitude, south first (as the original
    tool assigned sorted values: lowest value to the southernmost cell)."""
    lats = [float(np.mean(r[:-1, 1])) for r in rings]
    return list(np.argsort(lats, kind='stable'))


def inset_ring(ring, center, frac):
    """The bin ring shrunk toward its center by frac (0.08 = 8% of each vertex's
    distance), in the gnomonic plane: strictly inside the bin, so neighbouring bins'
    inset rings never touch."""
    lon0, lat0 = center
    xy = _to_plane(np.asarray(ring, dtype=float), lon0, lat0)
    return _to_lonlat(xy * (1.0 - frac), lon0, lat0)
