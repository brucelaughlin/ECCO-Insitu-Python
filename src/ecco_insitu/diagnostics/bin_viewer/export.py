"""
Data for the geodesic bin viewer, a web copy of the original ``geodesic_binning`` plotting app.

Two kinds of files, written to one folder that is published as is (GitHub Pages, or any static host):

``bins.bin`` (gzip) + ``bins.json``
    The geodesic bin polygons, built exactly as the original binning step builds them
    (``geodesic_binning/binning/geodesic_binning_utilities_binning.py``): the nearest-vertex bins of the
    10242-vertex geodesic grid; each bin's polygon is the convex hull of the points of a 0.1 deg lon/lat grid
    that fall in it; bins spanning more than 100 deg of longitude are dropped (the original's dateline
    "map_span_bug"); each polygon is clipped to the ocean with Natural Earth's 10 m land, dropped if less than
    2 % of it is ocean, and reduced to its largest piece. Lon/lat degrees. Also each bin's representative
    point (for the bin-mean labels). The polygons don't depend on the data: built once.

``bins_N.bin``, ``bins_S.bin``, ``land_N.json``, ``land_S.json`` (gzip)
    The polar maps: the bins and land of each cap, in polar stereographic map units (:func:`build_polar_bins`).

``land.json``, ``bathy.json`` (gzip)
    Natural Earth 10 m land (simplified to 0.02 deg) and isobaths (200, 1000 ... 5000 m), for drawing.

``<class>/profiles.bin`` + ``<class>/<V>_<level>.bin`` (gzip) + ``index.json``
    Per obs class: for every profile its lon, lat, year and bin (sorted by bin, then year), and for each
    variable and depth level the anomaly V - V_clim of the profiles that have a value there (a profile
    index list and float16 values). The browser computes bin statistics for any time window from these.
    Only values the chain kept (weight > 0) are written; the original app used every non-missing value.
    ``index.json`` also holds each class's color limits, computed as the original computes them (quantiles of
    all values / of the bins' std, over all depth levels) but over the whole record, so the colors stay
    fixed while the time window changes.

Layouts (little endian):
    bins.bin:      uint32 n_bins, then per bin: uint16 bin id, uint16 n_vertices, float32 rep lon, rep lat,
                   n_vertices x (float32 lon, float32 lat)
    profiles.bin:  uint32 n, then byte planes of: float32 lon[n], float32 lat[n], uint16 year[n],
                   uint16 bin[n] as differences (sorted by bin, then year)
    <V>_<k>.bin:   uint32 m, then byte planes of: uint32 profile index[m] (into profiles.bin, ascending) as
                   differences, float16 anomaly[m]
    ("byte planes": all first bytes of the column, then all second bytes, ...; lossless; see FORMAT)

Usage::

    python -m ecco_insitu.diagnostics.bin_viewer.export --run <run>/yearly --out <folder> [--sources PFL ...]
    python -m ecco_insitu.diagnostics.bin_viewer.export --out <folder> --page-only      # after editing viewer.html
    python -m ecco_insitu.diagnostics.bin_viewer.export --out <folder> --upload <bucket> --endpoint <url>

The folder also gets ``index.html`` (the page) and ``vendor/`` (deck.gl, d3-delaunay), so it is the whole site.
See docs/bin_viewer.md.
"""
import argparse
import gzip
import re
import json
import time
from pathlib import Path

import netCDF4
import numpy as np

import ecco_insitu.config as config_loader
import ecco_insitu.diagnostics.common as common

VARS = ('T', 'S')
ANGULAR_PRECISION = 0.1        # deg, the original's grid for the bin polygons
MIN_OCEAN_FRACTION = 0.02      # the original's min_ocean_fraction
LAND_RESOLUTION = '10m'
QF_FACE = (0.05, 0.95)         # the original's quantiles for the face (mean) colors: all values, all levels
QF_EDGE = (0.05, 0.25, 0.5, 0.75, 0.95, 0.99)   # ... and for the edge (std) colors: bin std, all levels


def _sph2cart(lon_deg, lat_deg):
    lon, lat = np.radians(lon_deg), np.radians(lat_deg)
    return np.stack([np.cos(lat) * np.cos(lon), np.cos(lat) * np.sin(lon), np.sin(lat)], axis=-1)


def _geodesic_tree():
    from scipy.spatial import cKDTree
    cfg = config_loader.load()
    ll = np.loadtxt(config_loader.get(cfg, 'paths.geodesic_bin_csv'), delimiter=',')
    return cKDTree(_sph2cart(ll[:, 0], ll[:, 1])), len(ll)


def build_bins(out):
    """Bin polygons as the original builds them; writes bins.bin and bins.json."""
    import cartopy.io.shapereader as shpreader
    from scipy.spatial import ConvexHull
    from shapely.geometry import Polygon
    from shapely.ops import unary_union
    from shapely.prepared import prep

    t0 = time.time()
    tree, nbins = _geodesic_tree()
    lons = np.arange(-180, 180, ANGULAR_PRECISION)
    lats = np.arange(-90, 90, ANGULAR_PRECISION)
    LON, LAT = np.meshgrid(lons, lats)
    _, idx = tree.query(_sph2cart(LON.ravel(), LAT.ravel()))      # 0-based bin of each grid point
    LON, LAT = LON.ravel(), LAT.ravel()
    order = np.argsort(idx, kind='stable')
    idx_s, lon_s, lat_s = idx[order], LON[order], LAT[order]
    starts = np.searchsorted(idx_s, np.arange(nbins + 1))
    print(f'  bin grid: {time.time() - t0:.0f}s')

    land = unary_union(list(shpreader.Reader(
        shpreader.natural_earth(resolution=LAND_RESOLUTION, category='physical', name='land')).geometries()))
    land_p = prep(land)
    print(f'  land polygon: {time.time() - t0:.0f}s')

    chunks, kept, dropped = [], 0, {'span': 0, 'land': 0, 'small': 0, 'tiny': 0}
    for b in range(nbins):
        lo, la = lon_s[starts[b]:starts[b + 1]], lat_s[starts[b]:starts[b + 1]]
        if len(lo) < 3:
            dropped['tiny'] += 1
            continue
        if lo.max() - lo.min() > 100:
            dropped['span'] += 1
            continue
        pts = np.c_[lo, la]
        try:
            hull = pts[ConvexHull(pts).vertices]
        except Exception:
            dropped['tiny'] += 1
            continue
        poly = Polygon(hull)
        if land_p.intersects(poly):
            ocean = poly.difference(land)
        else:
            ocean = poly
        if ocean.is_empty:
            dropped['land'] += 1
            continue
        if ocean.area / poly.area < MIN_OCEAN_FRACTION:
            dropped['small'] += 1
            continue
        if ocean.geom_type == 'MultiPolygon':
            ocean = max(ocean.geoms, key=lambda g: g.area)
        if ocean.geom_type != 'Polygon':
            dropped['small'] += 1
            continue
        xy = np.asarray(ocean.exterior.coords, np.float32)
        rp = ocean.representative_point()
        chunks += [np.array([b + 1, len(xy)], '<u2').tobytes(), np.array([rp.x, rp.y], '<f4').tobytes(),
                   xy.astype('<f4').tobytes()]
        kept += 1
    blob = np.uint32(kept).tobytes() + b''.join(chunks)
    (out / 'bins.bin').write_bytes(gzip.compress(blob, 9))
    (out / 'bins.json').write_text(json.dumps({
        'n_bins_grid': int(nbins), 'n_bins_kept': kept, 'dropped': dropped,
        'angular_precision': ANGULAR_PRECISION, 'min_ocean_fraction': MIN_OCEAN_FRACTION,
        'land': f'Natural Earth {LAND_RESOLUTION} land'}))
    print(f'  bins: {kept} kept, dropped {dropped}; {len(blob) / 1e6:.1f} MB raw; {time.time() - t0:.0f}s')


# ---- file layout (FORMAT 2) ----
# Lossless: every number is stored exactly; only the byte order changes, so gzip packs it better. Each
# column is split into byte planes (all first bytes, then all second bytes, ...); sorted integer columns
# (profile index, bin) are stored as differences from the previous value. About half the size of the
# plain layout (FORMAT 1: whole little-endian numbers).
FORMAT = 2


def _planes(a):
    """Little-endian byte planes of a 1-D array: plane 0 = lowest bytes of every value, ..."""
    return np.ascontiguousarray(a.astype(a.dtype.newbyteorder('<')).view(np.uint8).reshape(-1, a.itemsize).T).tobytes()


def _unplanes(buf, n, dtype):
    dt = np.dtype(dtype).newbyteorder('<')
    return np.frombuffer(buf, np.uint8, n * dt.itemsize).reshape(dt.itemsize, n).T.copy().view(dt).ravel()


def _delta(a):
    return np.diff(a.astype(np.int64), prepend=0)


def encode_profiles(lon, lat, year, bins):
    """profiles.bin: uint32 n; byte planes of f32 lon, f32 lat, u16 year, u16 delta(bin) (bins ascending)."""
    n = len(lon)
    return gzip.compress(np.uint32(n).tobytes() + _planes(lon.astype('<f4')) + _planes(lat.astype('<f4'))
                         + _planes(year.astype('<u2')) + _planes(_delta(bins).astype('<u2')), 9)


def decode_profiles(blob):
    """Inverse of :func:`encode_profiles`: lon, lat, year, bin arrays."""
    b = gzip.decompress(blob)
    n = int(np.frombuffer(b, '<u4', 1)[0])
    o, out = 4, []
    for dt in ('<f4', '<f4', '<u2', '<u2'):
        w = np.dtype(dt).itemsize * n
        out.append(_unplanes(b[o:o + w], n, dt))
        o += w
    out[3] = np.cumsum(out[3].astype(np.int64)).astype('<u2')
    return out


def encode_level(ix, v):
    """<V>_<k>.bin: uint32 m; byte planes of u32 delta(profile index) (ascending), f16 anomaly."""
    return gzip.compress(np.uint32(len(ix)).tobytes() + _planes(_delta(ix).astype('<u4'))
                         + _planes(v.astype('<f2')), 9)


def decode_level(blob):
    """Inverse of :func:`encode_level`: profile index and float16 anomaly arrays."""
    b = gzip.decompress(blob)
    m = int(np.frombuffer(b, '<u4', 1)[0])
    ix = np.cumsum(_unplanes(b[4:4 + 4 * m], m, '<u4').astype(np.int64)).astype('<u4')
    return ix, _unplanes(b[4 + 4 * m:4 + 6 * m], m, '<f2')


def convert_format1(out):
    """Rewrite a FORMAT 1 folder as FORMAT 2 in place, checking each file decodes to the same numbers."""
    out = Path(out)
    index = json.loads((out / 'index.json').read_text())
    if index.get('format', 1) >= 2:
        print('already format 2')
        return
    before = after = 0
    for src in index['classes']:
        for f in sorted((out / src).glob('*.bin')):
            old = f.read_bytes()
            raw = gzip.decompress(old)
            n = int(np.frombuffer(raw, '<u4', 1)[0])
            if f.name == 'profiles.bin':
                cols = [np.frombuffer(raw, dt, n, 4 + off) for dt, off in
                        (('<f4', 0), ('<f4', 4 * n), ('<u2', 8 * n), ('<u2', 10 * n))]
                new = encode_profiles(*cols)
                back = decode_profiles(new)
            else:
                cols = [np.frombuffer(raw, '<u4', n, 4), np.frombuffer(raw, '<f2', n, 4 + 4 * n)]
                new = encode_level(*cols)
                back = decode_level(new)
            if not all(a.tobytes() == c.tobytes() for a, c in zip(cols, back)):   # bit-for-bit, NaNs included
                raise RuntimeError(f'{f}: round trip differs; nothing more converted')
            f.write_bytes(new)
            before += len(old)
            after += len(new)
        print(f'  {src}: done')
    index['format'] = FORMAT
    (out / 'index.json').write_text(json.dumps(index, indent=1))
    print(f'format 2: {before / 1e6:.0f} MB -> {after / 1e6:.0f} MB, every file checked bit for bit')


def _bin_std(bins, v):
    """np.std (ddof 0) of v per bin, for the bins that have values."""
    u, inv, n = np.unique(bins, return_inverse=True, return_counts=True)
    s1, s2 = np.bincount(inv, v), np.bincount(inv, v * v)
    return np.sqrt(np.maximum(s2 / n - (s1 / n) ** 2, 0))


def _kept_bins(out):
    b = gzip.decompress((out / 'bins.bin').read_bytes())
    n, o, ids = int(np.frombuffer(b, '<u4', 1)[0]), 4, []
    for _ in range(n):
        bid, nv = np.frombuffer(b, '<u2', 2, o)
        ids.append(int(bid))
        o += 4 + 8 + 8 * int(nv)
    return np.array(ids)


# ---- polar maps ----
# Polar stereographic (true scale at the pole), on a sphere; units: km / KM_PER_UNIT, so one unit is about
# one degree at the equator and the page's zoom thresholds mean the same distance on every map.
EARTH_KM = 6371.0
KM_PER_UNIT = 111.195
POLAR_MIN_LAT = 20         # a cap's bins: those with any point poleward of this (the square map's corners
                           # reach about 23 deg when the 40 deg circle fills its width)


def polar_xy(lon, lat, hemi):
    """Polar stereographic x, y (map units) of lon/lat degrees; hemi 'N' or 'S' (as cartopy's Stereographic)."""
    lo, la = np.radians(lon), np.radians(lat)
    if hemi == 'N':
        r = 2 * EARTH_KM * np.tan(np.pi / 4 - la / 2)
        return r * np.sin(lo) / KM_PER_UNIT, -r * np.cos(lo) / KM_PER_UNIT
    r = 2 * EARTH_KM * np.tan(np.pi / 4 + la / 2)
    return r * np.sin(lo) / KM_PER_UNIT, r * np.cos(lo) / KM_PER_UNIT


def build_polar_bins(out):
    """The bins of both polar caps, built as the original builds bins but in polar stereographic space.

    Same steps as :func:`build_bins` (0.1 deg grid points of each bin, convex hull, clip to land, drop below
    2 % ocean, largest piece), with the grid points projected first, so bins around the poles and across the
    dateline are whole (the lon/lat map drops those). Writes bins_N.bin, bins_S.bin (same layout as
    bins.bin, in map units) and land_N.json, land_S.json (projected land rings, for drawing).
    """
    import cartopy.crs as ccrs
    import cartopy.io.shapereader as shpreader
    from scipy.spatial import ConvexHull
    import shapely
    from shapely import affinity
    from shapely.geometry import Polygon, box
    from shapely.ops import unary_union
    from shapely.prepared import prep

    t0 = time.time()
    tree, nbins = _geodesic_tree()
    lons = np.arange(-180, 180, ANGULAR_PRECISION)
    lats = np.arange(-90, 90 + ANGULAR_PRECISION / 2, ANGULAR_PRECISION)
    LON, LAT = np.meshgrid(lons, lats)
    LON, LAT = LON.ravel(), LAT.ravel()
    _, idx = tree.query(_sph2cart(LON, LAT))
    globe = ccrs.Globe(ellipse=None, semimajor_axis=EARTH_KM * 1000, semiminor_axis=EARTH_KM * 1000)
    land_ll = list(shpreader.Reader(shpreader.natural_earth(resolution=LAND_RESOLUTION, category='physical',
                                                            name='land')).geometries())
    for hemi, sign in (('N', 1), ('S', -1)):
        crs = ccrs.Stereographic(central_latitude=90 * sign, central_longitude=0, globe=globe)
        # land in map units: the part of the hemisphere poleward of 10 deg, projected by cartopy (handles the pole)
        cap = box(-180, 10, 180, 90) if sign > 0 else box(-180, -90, 180, -10)
        parts = [g.intersection(cap) for g in land_ll]
        proj = [crs.project_geometry(g, ccrs.PlateCarree()) for g in parts if not g.is_empty]
        land = unary_union([shapely.make_valid(g).buffer(0) for g in proj if not g.is_empty])
        land = affinity.scale(land, 1 / (KM_PER_UNIT * 1000), 1 / (KM_PER_UNIT * 1000), origin=(0, 0))
        land_p = prep(land)
        sel = (LAT * sign) >= POLAR_MIN_LAT - 2
        ids = np.unique(idx[sel])
        order = np.argsort(idx, kind='stable')
        idx_s = idx[order]
        starts = np.searchsorted(idx_s, np.arange(nbins + 1))
        chunks, kept = [], 0
        for b in ids:
            j = order[starts[b]:starts[b + 1]]
            if (LAT[j] * sign).max() < POLAR_MIN_LAT or len(j) < 3:
                continue
            x, y = polar_xy(LON[j], LAT[j], hemi)
            pts = np.c_[x, y]
            try:
                poly = Polygon(pts[ConvexHull(pts).vertices])
            except Exception:
                continue
            ocean = poly.difference(land) if land_p.intersects(poly) else poly
            if ocean.is_empty or ocean.area / poly.area < MIN_OCEAN_FRACTION:
                continue
            if ocean.geom_type == 'MultiPolygon':
                ocean = max(ocean.geoms, key=lambda g: g.area)
            if ocean.geom_type != 'Polygon':
                continue
            xy = np.asarray(ocean.exterior.coords, np.float32)
            rp = ocean.representative_point()
            chunks += [np.array([b + 1, len(xy)], '<u2').tobytes(), np.array([rp.x, rp.y], '<f4').tobytes(),
                       xy.astype('<f4').tobytes()]
            kept += 1
        (out / f'bins_{hemi}.bin').write_bytes(gzip.compress(np.uint32(kept).tobytes() + b''.join(chunks), 9))
        rings = []
        for g in (land.geoms if land.geom_type == 'MultiPolygon' else [land]):
            g = g.simplify(0.02)
            if g.is_empty:
                continue
            for q in (g.geoms if g.geom_type == 'MultiPolygon' else [g]):
                rings.append([np.round(np.asarray(r.coords), 3).tolist() for r in (q.exterior, *q.interiors)])
        (out / f'land_{hemi}.json').write_bytes(gzip.compress(json.dumps(rings, separators=(',', ':')).encode(), 9))
        print(f'  polar {hemi}: {kept} bins; {time.time() - t0:.0f}s')


def write_land(out):
    """Natural Earth 10 m land, simplified to 0.02 deg, as rings of lon/lat (land.json, gzip)."""
    import cartopy.io.shapereader as shpreader
    polys = []
    for g in shpreader.Reader(shpreader.natural_earth(resolution=LAND_RESOLUTION, category='physical',
                                                      name='land')).geometries():
        g = g.simplify(0.02)
        for p in (g.geoms if g.geom_type == 'MultiPolygon' else [g]):
            polys.append([np.round(np.asarray(r.coords), 3).tolist() for r in (p.exterior, *p.interiors)])
    (out / 'land.json').write_bytes(gzip.compress(json.dumps(polys, separators=(',', ':')).encode(), 9))


BATHY_LEVELS = {200: 'K_200', 1000: 'J_1000', 2000: 'I_2000', 3000: 'H_3000', 4000: 'G_4000', 5000: 'F_5000'}


def write_bathy(out):
    """Isobaths (Natural Earth 10 m bathymetry outlines, simplified to 0.03 deg) as lon/lat paths (bathy.json, gzip).

    The outlines of Natural Earth's "deeper than d" polygons, split at their seams (edges over 1 deg and the
    dateline cut), without rings smaller than 0.3 deg.
    """
    import cartopy.io.shapereader as shpreader
    levels = []
    for d, name in BATHY_LEVELS.items():
        paths = []
        for g in shpreader.Reader(shpreader.natural_earth(resolution=LAND_RESOLUTION, category='physical',
                                                          name=f'bathymetry_{name}')).geometries():
            g = g.simplify(0.03)
            for p in (g.geoms if g.geom_type == 'MultiPolygon' else [g]):
                for r in (p.exterior, *p.interiors):
                    c = np.round(np.asarray(r.coords), 2)
                    # split at the polygon seams (long straight edges, the dateline cut) ...
                    jump = np.r_[True, np.hypot(*np.diff(c, axis=0).T) > 1.0]
                    for piece in np.split(c, np.nonzero(jump)[0][1:]):
                        piece = piece[np.abs(piece[:, 0]) < 179.99]
                        # ... and drop rings smaller than 0.3 deg (seamount clutter)
                        if len(piece) > 1 and np.ptp(piece, axis=0).max() >= 0.3:
                            paths.append(piece.tolist())
        levels.append({'depth': d, 'paths': paths})
    (out / 'bathy.json').write_bytes(gzip.compress(json.dumps(levels, separators=(',', ':')).encode(), 9))


def export_class(src, files, out, tree, kept_bins, estim=None, ref_name='model', ref_only=False):
    """Profiles and anomalies of one obs class: V - V_clim always, and V - V_model if ``estim`` (an
    :class:`ecco_insitu.diagnostics.profile_analyses.compute.EstimIndex`) is given (in ``<class>/<ref_name>/``).
    ``ref_only``: write only that model's level files (a second model added to an existing folder; the profile
    order is recomputed the same way, so it matches the folder's profiles.bin)."""
    t0 = time.time()
    lon, lat, year, nlev = [], [], [], None
    refs = ((ref_name,) if ref_only else ('clim', ref_name)) if estim is not None else ('clim',)
    vals = {(r, V): [] for r in refs for V in VARS}   # per file: (n_prof, n_lev) float32 anomaly, NaN = none
    for f in files:
        with netCDF4.Dataset(f) as ds:
            ds.set_auto_mask(False)
            n = ds.dimensions['iPROF'].size
            if n == 0:
                continue
            lo = ds['prof_lon'][:].astype(np.float64)
            lo = np.where(lo > 180, lo - 360, lo)
            lon.append(lo); lat.append(ds['prof_lat'][:].astype(np.float64))
            year.append((ds['prof_YYYYMMDD'][:] // 10000).astype(np.int32))
            nlev = ds.dimensions['iDEPTH'].size
            model = estim.get(f, n, nlev) if estim is not None else None
            for V in VARS:
                if f'prof_{V}' not in ds.variables:
                    for r in refs:
                        vals[r, V].append(np.full((n, nlev), np.nan, np.float32))
                    continue
                v = ds[f'prof_{V}'][:].astype(np.float64)
                w = ds[f'prof_{V}weight'][:].astype(np.float64)
                ok = (w > 0) & (w < 1e20) & (np.abs(v) < 900)
                if 'clim' in refs:
                    c = ds[f'prof_{V}clim'][:].astype(np.float64)
                    vals['clim', V].append(np.where(ok & (np.abs(c) < 900), v - c, np.nan).astype(np.float32))
                if estim is not None:
                    # the model's value at the profile; 0 or missing = no model value (as the cost analyses)
                    m = model[V] if model is not None and model.get(V) is not None else np.full((n, nlev), np.nan)
                    good = ok & np.isfinite(m) & (m != 0) & (np.abs(m) < 900)
                    vals[ref_name, V].append(np.where(good, v - m, np.nan).astype(np.float32))
    lon, lat, year = np.concatenate(lon), np.concatenate(lat), np.concatenate(year)
    good = np.isfinite(lon) & np.isfinite(lat) & (np.abs(lat) <= 90)
    _, b = tree.query(_sph2cart(lon, lat))
    b = (b + 1).astype(np.int64)
    order = np.lexsort((year, b))
    order = order[good[order]]
    d = out / src
    d.mkdir(parents=True, exist_ok=True)
    n = len(order)
    if not ref_only:
        (d / 'profiles.bin').write_bytes(encode_profiles(lon[order], lat[order], year[order], b[order]))
    bo, kept = b[order], np.isin(b[order], kept_bins)
    info = {'n_profiles': int(n), 'years': [int(year[order].min()), int(year[order].max())]}
    for r in refs:
        dr = d if r == 'clim' else d / r       # climatology files at the top (as before), others in a subfolder
        dr.mkdir(exist_ok=True)
        levels, scales = {}, {}
        for V in VARS:
            A = np.concatenate(vals[r, V])[order]
            have, allv, stds = [], [], []
            for k in range(nlev):
                col = A[:, k]
                i = np.nonzero(np.isfinite(col))[0]
                if not len(i):
                    continue
                (dr / f'{V}_{k:02d}.bin').write_bytes(encode_level(i, col[i]))
                have.append(k)
                allv.append(col[i][kept[i]])
                stds.append(_bin_std(bo[i][kept[i]], col[i][kept[i]].astype(np.float64)))
            levels[V] = have
            if have:
                allv, stds = np.concatenate(allv), np.concatenate(stds)
                scales[V] = {'face': [float(x) for x in np.quantile(allv, QF_FACE)],
                             'edge': [float(x) for x in np.quantile(stds, QF_EDGE)]}
        if r == 'clim':
            info.update(levels=levels, scales=scales)
        else:
            info.setdefault('refs', {})[r] = {'levels': levels, 'scales': scales}
        print(f'  {src} vs {r}: levels T {len(levels["T"])} S {len(levels["S"])}')
    size = sum(f.stat().st_size for f in d.rglob('*') if f.is_file()) / 1e6
    print(f'  {src}: {n:,} profiles, {size:.0f} MB, {time.time() - t0:.0f}s')
    return info


def platform_of(src, descr):
    """The platform a profile comes from, from its ``prof_descr`` ('' if none):

    - PFL: the Argo float (WMO number), e.g. 'ARGO-PFL-aoml//D5902328-1' -> '5902328';
    - ITP: the ITP number, e.g. 'ITP2_itp31grd0001_nk' -> 'itp31';
    - MEOP: the tag deployment, e.g. 'MEOP-D00050144-4/hd2bb-9340' -> 'D00050144 (hd2bb)';
    - MRB_WOD: the mooring number, e.g. 'WOD-MRB//15116734_23005_RESEAR' -> '23005' (99999 = unknown);
    - other WOD classes: the platform / ship code after the cast number, e.g. 'WOD-XBT//12284221:ONFC_2120234'
      -> 'ONFC', 'WOD-GLD//10880467:CO_61859_200' -> 'CO' (coarse: a ship or program, not one instrument).
    """
    d = descr.strip()
    if src == 'PFL':
        m = re.search(r'//[DRA]?(\d+)-', d)
        return m.group(1) if m else ''
    if src == 'ITP':
        m = re.search(r'itp(\d+)', d)
        return f'itp{m.group(1)}' if m else ''
    if src == 'MEOP':
        m = re.match(r'MEOP-(D\d+)-\d+/([^-\s]+)', d)
        return f'{m.group(1)} ({m.group(2)})' if m else ''
    if src == 'MRB_WOD':
        m = re.search(r'//\d+_(\d+)_', d)
        return m.group(1) if m and m.group(1) != '99999' else ''
    m = re.search(r'//\d+:([^_]*)', d)
    return m.group(1).strip() if m and m.group(1).strip() else ''


def write_platforms(run_dir, out, sources=None):
    """``<class>/platforms.bin`` (+ ``platforms.json``): per profile, in profiles.bin's order, the chain's geodesic bin
    (``prof_bin_id_a``, the bins of the cost analyses) and the platform (:func:`platform_of`), for the cost explorer's
    platform tables and the viewer's platform highlight.

    Reads the yearly files in the same order as :func:`export_class` and checks each class's lon/lat against
    profiles.bin, so the two files line up profile for profile. Layout: uint32 n; byte planes of uint16 chain
    bin, uint32 platform index (0xFFFFFFFF = none); platforms.json: the platform names, then the label of the
    platform kind for this class.
    """
    run_dir, out = common.yearly_dir(run_dir), Path(out)
    tree, _ = _geodesic_tree()
    groups = {}
    for f in sorted(run_dir.rglob('*.nc')):
        src, _ = common.classify(f, run_dir)
        if (not sources or src in sources) and (out / src / 'profiles.bin').exists():
            groups.setdefault(src, []).append(f)
    kinds = {'PFL': 'Argo float (WMO)', 'ITP': 'ITP', 'MEOP': 'tag deployment', 'MRB_WOD': 'mooring',
             'XBT_WOD': 'ship / platform code', 'CTD_WOD': 'ship / platform code', 'GLD_WOD': 'glider / program code'}
    for src, files in sorted(groups.items()):
        lon, lat, year, cb, plat = [], [], [], [], []
        for f in files:
            with netCDF4.Dataset(f) as ds:
                ds.set_auto_mask(False)
                if ds.dimensions['iPROF'].size == 0:
                    continue
                lo = ds['prof_lon'][:].astype(np.float64)
                lon.append(np.where(lo > 180, lo - 360, lo)); lat.append(ds['prof_lat'][:].astype(np.float64))
                year.append((ds['prof_YYYYMMDD'][:] // 10000).astype(np.int32))
                cb.append(ds['prof_bin_id_a'][:].astype(np.int64))
                plat += [platform_of(src, b''.join(x).decode(errors='replace')) for x in ds['prof_descr'][:]]
        lon, lat, year, cb = map(np.concatenate, (lon, lat, year, cb))
        plat = np.array(plat, dtype=object)
        good = np.isfinite(lon) & np.isfinite(lat) & (np.abs(lat) <= 90)
        _, b = tree.query(_sph2cart(lon, lat))
        order = np.lexsort((year, b + 1))
        order = order[good[order]]
        plon, plat_ = decode_profiles((out / src / 'profiles.bin').read_bytes())[:2]
        if len(plon) != len(order) or not (np.array_equal(plon, lon[order].astype('<f4'))
                                           and np.array_equal(plat_, lat[order].astype('<f4'))):
            raise RuntimeError(f'{src}: profile order differs from profiles.bin; re-export the class first')
        names, idx = np.unique(plat[order].astype(str), return_inverse=True)
        idx = idx.astype(np.uint32)
        if '' in names:                          # 'no platform' -> 0xFFFFFFFF, and drop '' from the names
            e = int(np.nonzero(names == '')[0][0])
            idx = np.where(idx == e, 0xFFFFFFFF, np.where(idx > e, idx - 1, idx)).astype(np.uint32)
            names = names[names != '']
        n = len(order)
        (out / src / 'platforms.bin').write_bytes(gzip.compress(
            np.uint32(n).tobytes() + _planes(np.clip(cb[order], 0, 65535).astype('<u2')) + _planes(idx.astype('<u4')), 9))
        (out / src / 'platforms.json').write_bytes(gzip.compress(json.dumps(
            {'kind': kinds.get(src, 'platform'), 'names': names.tolist()}).encode(), 9))
        none = int((idx == 0xFFFFFFFF).sum())
        print(f'  {src}: {n:,} profiles, {len(names):,} platforms ({kinds.get(src, "platform")}); '
              f'{none:,} without one; chain bin = viewer bin for {np.mean(cb[order] == b[order] + 1):.0%}')


# ---- profile detail (for the viewer's profile panel) ----
DETAIL_CHUNK = 2000            # profiles per detail file (consecutive in profiles.bin's order: by bin, then year)
DETAIL_SCALE = {'obs': 1000, 'clim': 1000, 'model': 1000, 'sig': 10000}   # stored as integers: 0.001, sigma 0.0001
# per value, a mask byte: bit 0 climatology present, bit 1 model present, bit 2 sigma present


def _ordered_rows(src, files, tree):
    """(file index, row) of every profile in profiles.bin's order (as export_class orders them), and lon/lat."""
    lon, lat, year, fi, ri = [], [], [], [], []
    for i, f in enumerate(files):
        with netCDF4.Dataset(f) as ds:
            ds.set_auto_mask(False)
            n = ds.dimensions['iPROF'].size
            if n == 0:
                continue
            lo = ds['prof_lon'][:].astype(np.float64)
            lon.append(np.where(lo > 180, lo - 360, lo)); lat.append(ds['prof_lat'][:].astype(np.float64))
            year.append((ds['prof_YYYYMMDD'][:] // 10000).astype(np.int32))
            fi.append(np.full(n, i, np.int32)); ri.append(np.arange(n, dtype=np.int32))
    lon, lat, year, fi, ri = map(np.concatenate, (lon, lat, year, fi, ri))
    good = np.isfinite(lon) & np.isfinite(lat) & (np.abs(lat) <= 90)
    _, b = tree.query(_sph2cart(lon, lat))
    order = np.lexsort((year, b + 1))
    order = order[good[order]]
    return fi[order], ri[order], lon[order], lat[order]


def write_details(run_dir, out, sources=None, estim_dir=None):
    """``detail/<class>/c<NNNN>.bin``: everything about each profile, for the viewer's profile panel, in chunks of
    DETAIL_CHUNK profiles in profiles.bin's order (so profile j is in chunk j // DETAIL_CHUNK).

    Per chunk (gzip; little endian): uint32 P (profiles); uint32 date[P] (YYYYMMDD), uint32 time[P] (HHMMSS); then per
    variable T, S: uint8 n[P] (values per profile), and for the n·P values in profile order: uint8 level, uint8 mask,
    uint16 QC code (step 06's prof_Vweight_code; 0 = kept), and int32 obs, clim, model, sigma as integers
    (DETAIL_SCALE), each stored as differences from the previous value of the same profile, as byte planes; then the
    profiles' prof_descr, one per line (UTF-8). Values are those of the yearly files, rounded to 0.001 (sigma 0.0001).
    """
    from ecco_insitu.diagnostics.profile_analyses.compute import EstimIndex
    run_dir, out = common.yearly_dir(run_dir), Path(out)
    tree, _ = _geodesic_tree()
    estim = EstimIndex(estim_dir, run_dir) if estim_dir else None
    groups = {}
    for f in sorted(run_dir.rglob('*.nc')):
        src, _ = common.classify(f, run_dir)
        if (not sources or src in sources) and (out / src / 'profiles.bin').exists():
            groups.setdefault(src, []).append(f)
    info = {}
    for src, files in sorted(groups.items()):
        t0 = time.time()
        fi, ri, lon, lat = _ordered_rows(src, files, tree)
        plon, plat = decode_profiles((out / src / 'profiles.bin').read_bytes())[:2]
        if len(plon) != len(fi) or not (np.array_equal(plon, lon.astype('<f4')) and np.array_equal(plat, lat.astype('<f4'))):
            raise RuntimeError(f'{src}: profile order differs from profiles.bin; re-export the class first')
        # read every file once: per (file, row) the profile's values (kept as flat arrays per file, with offsets)
        per = []
        for i, f in enumerate(files):
            with netCDF4.Dataset(f) as ds:
                ds.set_auto_mask(False)
                n = ds.dimensions['iPROF'].size
                rec = {'date': ds['prof_YYYYMMDD'][:].astype(np.uint32) if n else np.zeros(0, np.uint32),
                       'time': ds['prof_HHMMSS'][:].astype(np.uint32) if n and 'prof_HHMMSS' in ds.variables else np.zeros(n, np.uint32),
                       'descr': [b''.join(x).decode(errors='replace').strip() for x in ds['prof_descr'][:]] if n else []}
                m = estim.get(f, n, ds.dimensions['iDEPTH'].size) if (estim is not None and n) else None
                for V in VARS:
                    if not n or f'prof_{V}' not in ds.variables:
                        rec[V] = None
                        continue
                    v = ds[f'prof_{V}'][:].astype(np.float64)
                    have = np.isfinite(v) & (np.abs(v) < 900)
                    def col(name, scale):
                        if name not in ds.variables:
                            return np.zeros(int(have.sum()), np.int32), np.zeros(int(have.sum()), bool)
                        a = ds[name][:].astype(np.float64)[have]
                        ok = np.isfinite(a) & (np.abs(a) < 900)
                        return np.where(ok, np.round(np.nan_to_num(a) * scale), 0).astype(np.int32), ok
                    obs = np.round(v[have] * DETAIL_SCALE['obs']).astype(np.int32)
                    clim, cok = col(f'prof_{V}clim', DETAIL_SCALE['clim'])
                    sig, sok = col(f'prof_{V}uncertainty', DETAIL_SCALE['sig'])
                    if m is not None and m.get(V) is not None:
                        mv = np.asarray(m[V], np.float64)[have]
                        mok = np.isfinite(mv) & (mv != 0) & (np.abs(mv) < 900)
                        mod = np.where(mok, np.round(np.nan_to_num(mv) * DETAIL_SCALE['model']), 0).astype(np.int32)
                    else:
                        mod, mok = np.zeros(len(obs), np.int32), np.zeros(len(obs), bool)
                    code = ds[f'prof_{V}weight_code'][:][have].astype(np.uint16) if f'prof_{V}weight_code' in ds.variables \
                        else np.zeros(len(obs), np.uint16)
                    cnt = have.sum(1)
                    rec[V] = {'cnt': cnt, 'off': np.r_[0, np.cumsum(cnt)], 'lev': np.nonzero(have)[1].astype(np.uint8),
                              'mask': (cok * 1 + mok * 2 + sok * 4).astype(np.uint8), 'code': code,
                              'obs': obs, 'clim': clim, 'model': mod, 'sig': sig}
                per.append(rec)
        d = out / 'detail' / src
        d.mkdir(parents=True, exist_ok=True)
        n_all, nch, size = len(fi), 0, 0
        for c0 in range(0, n_all, DETAIL_CHUNK):
            sl = slice(c0, min(c0 + DETAIL_CHUNK, n_all))
            F, Rw = fi[sl], ri[sl]
            P = len(F)
            date = np.array([per[f]['date'][r] for f, r in zip(F, Rw)], '<u4')
            tim = np.array([per[f]['time'][r] for f, r in zip(F, Rw)], '<u4')
            blob = [np.uint32(P).tobytes(), _planes(date), _planes(tim)]
            for V in VARS:
                cnts, cols = np.zeros(P, np.uint8), {k: [] for k in ('lev', 'mask', 'code', 'obs', 'clim', 'model', 'sig')}
                for q, (f, r) in enumerate(zip(F, Rw)):
                    rv = per[f][V]
                    if rv is None:
                        continue
                    a, b2 = rv['off'][r], rv['off'][r + 1]
                    cnts[q] = b2 - a
                    for k in cols:
                        cols[k].append(rv[k][a:b2])
                cat = {k: (np.concatenate(v) if v else np.zeros(0, np.int32 if k in DETAIL_SCALE else (np.uint16 if k == 'code' else np.uint8)))
                       for k, v in cols.items()}
                first = np.r_[0, np.cumsum(cnts.astype(np.int64))[:-1]][cnts > 0]
                def delta(x):
                    dd = np.diff(x.astype(np.int64), prepend=0)
                    dd[first] = x[first]
                    return dd.astype('<i4')
                blob += [cnts.tobytes(), cat['lev'].astype('u1').tobytes(), cat['mask'].astype('u1').tobytes(),
                         _planes(cat['code'].astype('<u2'))] + [_planes(delta(cat[k])) for k in ('obs', 'clim', 'model', 'sig')]
            blob.append('\n'.join(per[f]['descr'][r] for f, r in zip(F, Rw)).encode())
            data = gzip.compress(b''.join(blob), 9)
            (d / f'c{c0 // DETAIL_CHUNK:04d}.bin').write_bytes(data)
            nch += 1
            size += len(data)
        info[src] = {'chunk': DETAIL_CHUNK, 'n_chunks': nch, 'bytes': size}
        print(f'  {src}: {nch} detail chunks, {size / 1e6:.0f} MB, {time.time() - t0:.0f}s')
        del per
    ip = out / 'index.json'
    index = json.loads(ip.read_text())
    index.setdefault('detail', {}).update({k: {**v, 'scale': DETAIL_SCALE} for k, v in info.items()})
    ip.write_text(json.dumps(index, indent=1))


def decode_detail(blob):
    """Inverse of one detail chunk: dict with date, time, descr and per variable the values per profile (floats)."""
    b = gzip.decompress(blob)
    P = int(np.frombuffer(b, '<u4', 1)[0])
    o = 4
    out = {'date': _unplanes(b[o:o + 4 * P], P, '<u4')}; o += 4 * P
    out['time'] = _unplanes(b[o:o + 4 * P], P, '<u4'); o += 4 * P
    for V in VARS:
        cnt = np.frombuffer(b, 'u1', P, o).astype(np.int64); o += P
        n = int(cnt.sum())
        lev = np.frombuffer(b, 'u1', n, o); o += n
        mask = np.frombuffer(b, 'u1', n, o); o += n
        code = _unplanes(b[o:o + 2 * n], n, '<u2'); o += 2 * n
        first = np.r_[0, np.cumsum(cnt)[:-1]]
        vals = {}
        for k in ('obs', 'clim', 'model', 'sig'):
            d = _unplanes(b[o:o + 4 * n], n, '<i4').astype(np.int64); o += 4 * n
            # undo the per-profile differences: cumulative sum restarted at each profile's first value
            cs = np.cumsum(d)
            start = np.repeat(np.r_[0, cs][first[cnt > 0]] - 0, cnt[cnt > 0]) if n else np.zeros(0, np.int64)
            base = np.repeat(cs[first[cnt > 0]] - d[first[cnt > 0]], cnt[cnt > 0]) if n else np.zeros(0, np.int64)
            vals[k] = (cs - base) / DETAIL_SCALE[k]
        out[V] = {'cnt': cnt, 'lev': lev, 'mask': mask, 'code': code, **vals}
    out['descr'] = b[o:].decode().split('\n')
    return out


GUIDE_HTML = """
<p>A web copy of the original <code>geodesic_binning</code> plotting app. Everything is drawn as the original draws it
(colors, color limits, edge widths, the zoom threshold, the cells, dots and labels); the differences are listed below.</p>
<table>
<tr><th>Element</th><th>Meaning</th></tr>
<tr><td>bin face color</td><td>mean anomaly (V &minus; V<sub>clim</sub>) of the profiles in the bin, at this depth level and in the
  time window; PRGn, cut so white is 0; limits = 5 % / 95 % quantiles of all values of the variable (all depth levels)</td></tr>
<tr><td>bin face opacity</td><td>sqrt(count / largest count), between 0.4 and 1</td></tr>
<tr><td>bin edge color</td><td>standard deviation (ddof 0) of the anomalies in the bin; cividis, 25 % to 95 % quantiles of the bins' std</td></tr>
<tr><td>bin edge width</td><td>0.1 px &times; sqrt(profile count) (none for one profile), growing with the square root
  of the zoom factor; at most 25 % of the bin's own size (the original: 0.01 pt per profile, capped by the smallest bin;
  changed 2026-11-02 so busy regions, e.g. XBT shipping lanes, stand out zoomed out)</td></tr>
<tr><td>zoom level</td><td>sqrt(global area / view area); above 10 each bin with at most <i>max cells/bin</i> profiles splits into
  one cell per profile (Voronoi, 2 Lloyd steps), values sorted low to high from south to north; red dots are the true
  profile locations; numbers are the bin means</td></tr>
<tr><td>time</td><td>new: any window of whole years; "all years" is the whole record of the class</td></tr>
</table>
<p><b>Differences from the original.</b> (1) Only values the processing chain kept (weight &gt; 0) are shown; the original
used every non-missing value. (2) The color limits are computed per class over the whole record, so they stay fixed
while the time window changes. (3) The original showed one yearly file; this page shows any range of years, of all
files of the class. (4) If more than CELL_BUDGET cells would be in view, bins stay whole until you zoom further
(the browser's limit, not the original's). (5) Plain lon/lat map (as the original), drawn with WebGL (deck.gl).
(6) Faint isobaths (200, 1000, 2000, 3000, 4000, 5000 m; Natural Earth), not in the original; the 200 m line (the
shelf edge) is a little darker; each is labeled with its depth in meters along the line. The "isobaths" button
turns them off. (7) The view stays on the map: you can't pan past its edges or zoom out past it; "home" (or 9)
returns to the whole map. (8) Below the world map, an Arctic and a Southern Ocean map (polar stereographic, the 40&deg;
circle at the edges), with their own bins built the same way in polar coordinates, so bins around the poles and across
the dateline are whole; each zooms and pans on its own, with the same zoom threshold. (9) No keyboard shortcuts:
every control is a button. (10) Depth ranges as in the cost map explorer (standard bands, 100 m bands, 15 or 30
groups of levels), besides single levels: in a range each profile counts once, with the mean of its values at the
range's levels. (11) With model values exported (<code>--model-equivalents</code>), the anomaly can be relative to the
model (V &minus; V<sub>model</sub>) instead of the climatology; colors then use the model anomaly's own limits.
(12) A platform box (Argo float WMO number, ITP, tag, mooring, or WOD ship code) rings that platform's profiles in cyan at any
zoom; the cost map explorer's platform table links here with it filled in. (13) Zoomed in, <b>click a profile</b> (a red dot
or a ring) for its panel: the full profile (observation, climatology, model, neighbors in gray; obs &minus; reference with &plusmn;1&sigma;
and &plusmn;2&sigma;; levels the chain rejected marked &times;), its platform over time, and the bin's largest misfits. Values there are
the processed files' to 0.001.</p>
<p>Data: profiles.bin + one file per variable and depth level per class, written by
<code>ecco_insitu.diagnostics.bin_viewer.export</code>; loaded on demand and cached for the session.</p>
"""


def write_page(out):
    """index.html (viewer.html with the colormaps and guide filled in) and the vendored libraries."""
    import shutil
    from matplotlib import colormaps
    here = Path(__file__).parent
    lut = lambda name: ','.join(str(int(round(255 * c))) for i in range(256) for c in colormaps[name](i)[:3])
    html = (here / 'viewer.html').read_text()
    html = html.replace('__RDBU__', lut('RdBu_r')).replace('__CIVIDIS__', lut('cividis'))
    html = html.replace('__GUIDE__', GUIDE_HTML.replace('CELL_BUDGET', '200,000'))
    (out / 'index.html').write_text(html)
    shutil.copytree(here / 'vendor', out / 'vendor', dirs_exist_ok=True)


def add_model(run_dir, out, name, estim_dir, site_out, url=None, sources=None):
    """Add a second (third, ...) model to an existing viewer folder ``out``: its anomaly level files go to
    ``site_out/<out name>/<class>/<name>/`` (published as a site of their own), and ``out/index.json`` gets the
    model's levels, color scales, label and (``url``) where its files are. The page then offers it next to the
    climatology and the first model."""
    from ecco_insitu.diagnostics.profile_analyses.compute import EstimIndex
    run_dir, out, site_out = common.yearly_dir(run_dir), Path(out), Path(site_out)
    index = json.loads((out / 'index.json').read_text())
    groups = {}
    for f in sorted(run_dir.rglob('*.nc')):
        src, _ = common.classify(f, run_dir)
        if (not sources or src in sources) and src in index['classes']:
            groups.setdefault(src, []).append(f)
    tree, _ = _geodesic_tree()
    estim = EstimIndex(estim_dir, run_dir)
    dest = site_out / out.name
    dest.mkdir(parents=True, exist_ok=True)
    (dest / '.nojekyll').write_text('')
    for src, files in sorted(groups.items()):
        info = export_class(src, files, dest, tree, _kept_bins(out), estim, ref_name=name, ref_only=True)
        index['classes'][src].setdefault('refs', {})[name] = info['refs'][name]
    index.setdefault('refs', {'clim': 'climatology'})[name] = name
    if url:
        index.setdefault('ref_urls', {})[name] = url.rstrip('/')
    (out / 'index.json').write_text(json.dumps(index, indent=1))
    size = sum(p.stat().st_size for p in dest.rglob('*') if p.is_file()) / 1e6
    print(f'{name}: level files in {dest} ({size:.0f} MB); index.json updated')
    return dest


def run(run_dir, out, sources=None, rebuild_bins=False, estim_dir=None):
    """Write the whole viewer folder for one processed run: bins, land, isobaths, every class, the page."""
    run_dir = common.yearly_dir(run_dir)
    out = Path(out)
    out.mkdir(parents=True, exist_ok=True)
    if rebuild_bins or not (out / 'bins.bin').exists():
        print('Bins (once):')
        build_bins(out)
    if not (out / 'land.json').exists():
        write_land(out)
    if not (out / 'bathy.json').exists():
        write_bathy(out)
    if not (out / 'bins_N.bin').exists():
        build_polar_bins(out)
    kept_bins = _kept_bins(out)
    groups = {}
    for f in sorted(run_dir.rglob('*.nc')):
        src, _ = common.classify(f, run_dir)
        if not sources or src in sources:
            groups.setdefault(src, []).append(f)
    tree, _ = _geodesic_tree()
    estim = None
    if estim_dir:
        from ecco_insitu.diagnostics.profile_analyses.compute import EstimIndex
        estim = EstimIndex(estim_dir, run_dir)
    depth = None
    with netCDF4.Dataset(next(iter(groups.values()))[0]) as ds:
        depth = [float(x) for x in ds['prof_depth'][:]]
    index_path = out / 'index.json'
    index = json.loads(index_path.read_text()) if index_path.exists() else {'classes': {}}   # keeps report_url
    index.update({'format': FORMAT, 'run': str(run_dir), 'depth': depth, 'vars': list(VARS),
                  'units': {'T': '°C', 'S': 'psu'},
                  'created': time.strftime('%Y-%m-%dT%H:%M:%S')})
    print('Classes:')
    for src, files in sorted(groups.items()):
        index['classes'][src] = export_class(src, files, out, tree, kept_bins, estim)
    if estim_dir:
        name = Path(estim_dir).name if Path(estim_dir).name != 'netCDF' else Path(estim_dir).parent.name
        index['refs'] = {'clim': 'climatology', 'model': f'model ({name})'}
    else:
        index['refs'] = {'clim': 'climatology'}
    index_path.write_text(json.dumps(index, indent=1))
    total = sum(p.stat().st_size for p in out.rglob('*') if p.is_file()) / 1e6
    write_page(out)
    print(f'Done: {out} ({total:.0f} MB)')
    return out


GITHUB_LIMIT_MB = 1000         # GitHub Pages site limit (about 1 GB); every file must also be < 100 MB


def split_for_github(out, remote_name, url):
    """Move the model level files (``<class>/model/``) of ``out`` to ``<out>_<name>_site/<out name>/``, for a second GitHub Pages
    site, and record its address in index.json (``ref_urls``), so the page reads them from there.

    The two folders are then published separately (``--github`` on each, to two repositories), keeping each site
    under the 1 GB limit. Returns the new folder.
    """
    import shutil
    out = Path(out)
    dest = out.parent / f'{out.name}_{remote_name}_site'
    for d in sorted(out.glob('*/model')):
        target = dest / out.name / d.parent.name / 'model'
        target.parent.mkdir(parents=True, exist_ok=True)
        if target.exists():
            shutil.rmtree(target)
        shutil.move(str(d), str(target))
    (dest / out.name / '.nojekyll').write_text('')
    ip = out / 'index.json'
    index = json.loads(ip.read_text())
    index.setdefault('ref_urls', {})['model'] = url.rstrip('/')
    ip.write_text(json.dumps(index, indent=1))
    print(f'model files moved to {dest / out.name}; index.json points the page at {url}')
    return dest / out.name


def publish_github(out, remote, dry_run=False, landing=None, skip_detail=True):
    """Replace the viewer repository's main branch with this folder (one commit, force-pushed, no history).

    The repository (e.g. ``git@github.com:<user>/ecco-bin-viewer.git``) serves GitHub Pages from ``main``, root
    folder; the viewer is then at ``https://<user>.github.io/<repo>/<folder name>/``. One run at a time: the
    branch is replaced, so older runs go (they stay in the local folders).
    """
    import shutil
    import subprocess
    import tempfile
    out = Path(out).resolve()
    files = [f for f in out.rglob('*') if f.is_file() and not (skip_detail and f.relative_to(out).parts[:1] == ('detail',))]
    big = [f for f in files if f.stat().st_size >= 95e6]
    size = sum(f.stat().st_size for f in files) / 1e6
    if big:
        raise SystemExit(f'files over GitHub\'s 100 MB limit: {big[:3]}')
    if size > GITHUB_LIMIT_MB:
        raise SystemExit(f'{size:.0f} MB is over the GitHub Pages limit (~{GITHUB_LIMIT_MB} MB)')
    git = lambda *a, cwd: subprocess.run(['git', *a], cwd=cwd, check=True, capture_output=True, text=True)
    me = Path(__file__).resolve().parents[4]
    with tempfile.TemporaryDirectory() as tmp:
        work = Path(tmp) / 'site'
        work.mkdir()
        # the main site leaves out detail/ (published to its own sites by publish_details, which passes skip_detail=False)
        shutil.copytree(out, work / out.name, ignore=(lambda d, names: ['detail'] if skip_detail and Path(d) == out else []))
        (work / '.nojekyll').write_text('')
        if not (work / out.name / 'index.html').exists():
            target = landing or f'{out.name}/'
        else:
            target = f'{out.name}/'
        (work / 'index.html').write_text(
            f'<!DOCTYPE html><meta charset="utf-8"><meta http-equiv="refresh" content="0; url={target}">'
            f'<a href="{target}">Geodesic bin viewer ({out.name})</a>\n')
        (work / 'README.md').write_text(
            'Geodesic bin viewer (the in situ data viewer of the ECCO in situ profile reports). Made by '
            '`ecco-insitu bin-viewer --github`; this branch is replaced on every publish, so it has no history. '
            f'Run: {out.name}. Code: https://github.com/brucelaughlin/ECCO-Insitu-Python '
            '(src/ecco_insitu/diagnostics/bin_viewer).\n')
        print(f'{remote}: {out.name} ({size:.0f} MB)')
        if dry_run:
            print('dry run: nothing pushed')
            return
        git('init', '-q', '-b', 'main', cwd=work)
        git('add', '-A', cwd=work)
        name = git('config', 'user.name', cwd=me).stdout.strip()
        email = git('config', 'user.email', cwd=me).stdout.strip()
        git('-c', f'user.name={name}', '-c', f'user.email={email}', 'commit', '-q', '-m',
            f'Publish bin viewer: {out.name}', cwd=work)
        git('push', '-q', '--force', remote, 'main:main', cwd=work)
    print('pushed. First time only: repository -> Settings -> Pages -> Deploy from a branch -> main, / (root).')


def publish_details(out, remotes, dry_run=False):
    """Publish ``<out>/detail/<class>/`` to GitHub Pages sites: ``remotes`` maps a repository to the classes it holds
    (each site under ~1 GB). A class too big for one site is split by chunk number: write it as 'PFL:0-699' on one
    site and 'PFL:700-' on another. Writes the addresses into index.json (``detail_urls``: per class a list of
    [first chunk, last chunk, address]), so the page reads each chunk from its site; run before the main site."""
    import shutil
    import tempfile
    out = Path(out).resolve()
    ip = out / 'index.json'
    index = json.loads(ip.read_text())
    urls = {}
    for remote, classes in remotes.items():
        user, repo = re.match(r'.*[:/]([^/]+)/([^/]+?)(\.git)?$', remote).group(1, 2)
        url = f'https://{user}.github.io/{repo}/{out.name}'
        with tempfile.TemporaryDirectory() as tmp:
            site = Path(tmp) / out.name
            for spec in classes:
                c, _, rng = spec.partition(':')
                n = len(list((out / 'detail' / c).glob('c*.bin')))
                a, _, b = rng.partition('-') if rng else ('0', '', '')
                a, b = int(a or 0), int(b) if b else n - 1
                (site / 'detail' / c).mkdir(parents=True, exist_ok=True)
                for k in range(a, b + 1):
                    shutil.copy2(out / 'detail' / c / f'c{k:04d}.bin', site / 'detail' / c / f'c{k:04d}.bin')
                urls.setdefault(c, []).append([a, b, url])
            publish_github(site, remote, dry_run, landing=f'https://{user}.github.io/ecco-bin-viewer/{out.name}/', skip_detail=False)
    index['detail_urls'] = urls
    ip.write_text(json.dumps(index, indent=1))
    print('detail_urls:', urls)


def cli(argv=None):
    """``ecco-insitu bin-viewer``: build, convert, or publish (GitHub Pages) the viewer."""
    p = argparse.ArgumentParser(description='Write the geodesic bin viewer data.')
    p.add_argument('--add-model', nargs=3, metavar=('NAME', 'DIR', 'SITE_DIR'),
                   help='add a model to the existing --out folder: its name (e.g. V5r1), its model-equivalent folder, and '
                        'the folder for its own site (published with --github ... on that folder); needs --run')
    p.add_argument('--add-model-url', metavar='URL', help='with --add-model: the address its site will have')
    p.add_argument('--model-equivalents', metavar='DIR', default=None,
                   help='also V - model (the folder profile-analyze --model-equivalents reads)')
    p.add_argument('--run', help='processed run folder (or its yearly/)')
    p.add_argument('--out', required=True, help='output folder (published as is)')
    p.add_argument('--sources', nargs='*', default=None, help='obs classes (default: all)')
    p.add_argument('--rebuild-bins', action='store_true', help='rebuild bins.bin even if it exists')
    p.add_argument('--page-only', action='store_true', help='only rewrite index.html (after editing viewer.html)')
    p.add_argument('--bathy', action='store_true', help='(re)write bathy.json only')
    p.add_argument('--platforms', action='store_true', help='(re)write <class>/platforms.bin only (needs --run)')
    p.add_argument('--details', action='store_true', help='(re)write detail/<class>/ only (needs --run; '
                   'with --model-equivalents, also the model values)')
    p.add_argument('--polar', action='store_true', help='(re)build the polar maps\' bins and land only')
    p.add_argument('--report', metavar='URL', help='the profile report this viewer belongs to (its "← report" '
                   'button); kept in index.json')
    p.add_argument('--convert', action='store_true', help='rewrite an older (format 1) folder as format 2')
    p.add_argument('--github', metavar='REMOTE', help='then publish the folder to this GitHub Pages repository '
                   '(e.g. git@github.com:brucelaughlin/ecco-bin-viewer.git); replaces what is there')
    p.add_argument('--github-model', metavar='REMOTE', help='with --github: also publish <out>_model_site/ (made by '
                   '--split-model) to this repository (e.g. git@github.com:brucelaughlin/ecco-bin-viewer-model.git)')
    p.add_argument('--split-model', nargs=2, metavar=('NAME', 'URL'),
                   help='move the model files to <out>_NAME_site/, for a second GitHub Pages site at URL')
    p.add_argument('--github-detail', nargs='+', metavar='REMOTE=CLASS[:a-b],...', help='publish detail/ to these sites, '
                   'e.g. ...-detail-a.git=PFL:0-699 ...-detail-b.git=PFL:700- ...-detail-c.git=CTD_WOD,GLD_WOD,ITP,MEOP,MRB_WOD,XBT_WOD')
    p.add_argument('--dry-run', action='store_true', help='show what the publish would do')
    a = p.parse_args(argv)
    if a.convert:
        convert_format1(a.out)
    if a.add_model:
        if not a.run:
            p.error('--add-model needs --run')
        add_model(a.run, a.out, a.add_model[0], a.add_model[1], a.add_model[2], a.add_model_url, a.sources)
        return 0
    if a.split_model:
        split_for_github(a.out, *a.split_model)
    if a.bathy:
        write_bathy(Path(a.out))
    if a.platforms:
        if not a.run:
            p.error('--platforms needs --run')
        write_platforms(a.run, a.out, a.sources)
    if a.details:
        if not a.run:
            p.error('--details needs --run')
        write_details(a.run, a.out, a.sources, a.model_equivalents)
    if a.polar:
        build_polar_bins(Path(a.out))
    if a.report:
        ip = Path(a.out) / 'index.json'
        index = json.loads(ip.read_text())
        index['report_url'] = a.report
        ip.write_text(json.dumps(index, indent=1))
    if a.page_only or a.platforms or a.details or a.github_detail or a.bathy or a.polar or a.report or a.convert or a.split_model or (a.github and not a.run):
        write_page(Path(a.out))
    else:
        run(a.run, a.out, a.sources, a.rebuild_bins, a.model_equivalents)
    if a.github and not (a.page_only or a.platforms or a.bathy or a.polar or a.report or a.convert or a.split_model) and a.run:
        pass
    if a.github_detail:
        publish_details(a.out, {r.split('=')[0]: r.split('=')[1].split(',') for r in a.github_detail}, a.dry_run)
    if a.github:
        publish_github(Path(a.out), a.github, a.dry_run)
    if a.github_model:
        out = Path(a.out).resolve()
        model = out.parent / f'{out.name}_model_site' / out.name
        if not model.is_dir():
            p.error(f'{model} not found: run --split-model first')
        user, repo = re.match(r'.*[:/]([^/]+)/([^/]+?)(\.git)?$', a.github).group(1, 2)
        publish_github(model, a.github_model, a.dry_run, landing=f'https://{user}.github.io/{repo}/{out.name}/')


if __name__ == '__main__':
    cli()
