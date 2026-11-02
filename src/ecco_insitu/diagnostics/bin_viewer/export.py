"""
Data for the geodesic bin viewer, a web copy of the original ``geodesic_binning`` plotting app.

Two kinds of files, written to one folder that is uploaded as is (Cloudflare R2, or any static host):

``bins.bin`` (gzip) + ``bins.json``
    The geodesic bin polygons, built exactly as the original binning step builds them
    (``geodesic_binning/binning/geodesic_binning_utilities_binning.py``): the nearest-vertex bins of the
    10242-vertex geodesic grid; each bin's polygon is the convex hull of the points of a 0.1 deg lon/lat grid
    that fall in it; bins spanning more than 100 deg of longitude are dropped (the original's dateline
    "map_span_bug"); each polygon is clipped to the ocean with Natural Earth's 10 m land, dropped if less than
    2 % of it is ocean, and reduced to its largest piece. Lon/lat degrees. Also each bin's representative
    point (for the bin-mean labels). The polygons don't depend on the data: built once.

``land.json`` (gzip)
    Natural Earth 10 m land (simplified to 0.02 deg) for drawing.

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
    profiles.bin:  uint32 n, then columns: float32 lon[n], float32 lat[n], uint16 year[n], uint16 bin[n]
    <V>_<k>.bin:   uint32 m, then uint32 profile index[m] (into profiles.bin, ascending), float16 anomaly[m]

Usage::

    python -m ecco_insitu.diagnostics.bin_viewer.export --run <run>/yearly --out <folder> [--sources PFL ...]
    python -m ecco_insitu.diagnostics.bin_viewer.export --out <folder> --page-only      # after editing viewer.html
    python -m ecco_insitu.diagnostics.bin_viewer.export --out <folder> --upload <bucket> --endpoint <url>

The folder also gets ``index.html`` (the page) and ``vendor/`` (deck.gl, d3-delaunay), so it is the whole site.
See docs/bin_viewer.md.
"""
import argparse
import gzip
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


def export_class(src, files, out, tree, kept_bins):
    """Profiles and anomalies of one obs class."""
    t0 = time.time()
    lon, lat, year, nlev = [], [], [], None
    vals = {V: [] for V in VARS}            # per file: (n_prof, n_lev) float32 anomaly, NaN = none
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
            for V in VARS:
                if f'prof_{V}' not in ds.variables:
                    vals[V].append(np.full((n, nlev), np.nan, np.float32))
                    continue
                v = ds[f'prof_{V}'][:].astype(np.float64)
                c = ds[f'prof_{V}clim'][:].astype(np.float64)
                w = ds[f'prof_{V}weight'][:].astype(np.float64)
                ok = (w > 0) & (w < 1e20) & (np.abs(v) < 900) & (np.abs(c) < 900)
                vals[V].append(np.where(ok, v - c, np.nan).astype(np.float32))
    lon, lat, year = np.concatenate(lon), np.concatenate(lat), np.concatenate(year)
    good = np.isfinite(lon) & np.isfinite(lat) & (np.abs(lat) <= 90)
    _, b = tree.query(_sph2cart(lon, lat))
    b = (b + 1).astype(np.int64)
    order = np.lexsort((year, b))
    order = order[good[order]]
    d = out / src
    d.mkdir(parents=True, exist_ok=True)
    n = len(order)
    (d / 'profiles.bin').write_bytes(gzip.compress(
        np.uint32(n).tobytes() + lon[order].astype('<f4').tobytes() + lat[order].astype('<f4').tobytes()
        + year[order].astype('<u2').tobytes() + b[order].astype('<u2').tobytes(), 6))
    levels, scales = {}, {}
    bo, kept = b[order], np.isin(b[order], kept_bins)
    for V in VARS:
        A = np.concatenate(vals[V])[order]
        have, allv, stds = [], [], []
        for k in range(nlev):
            col = A[:, k]
            i = np.nonzero(np.isfinite(col))[0]
            if not len(i):
                continue
            (d / f'{V}_{k:02d}.bin').write_bytes(gzip.compress(
                np.uint32(len(i)).tobytes() + i.astype('<u4').tobytes() + col[i].astype('<f2').tobytes(), 6))
            have.append(k)
            allv.append(col[i][kept[i]])
            stds.append(_bin_std(bo[i][kept[i]], col[i][kept[i]].astype(np.float64)))
        levels[V] = have
        if have:
            allv, stds = np.concatenate(allv), np.concatenate(stds)
            scales[V] = {'face': [float(x) for x in np.quantile(allv, QF_FACE)],
                         'edge': [float(x) for x in np.quantile(stds, QF_EDGE)]}
    size = sum(p.stat().st_size for p in d.iterdir()) / 1e6
    print(f'  {src}: {n:,} profiles, levels T {len(levels["T"])} S {len(levels["S"])}, {size:.0f} MB, '
          f'{time.time() - t0:.0f}s')
    return {'n_profiles': int(n), 'levels': levels, 'scales': scales,
            'years': [int(year[order].min()), int(year[order].max())]}


GUIDE_HTML = """
<p>A web copy of the original <code>geodesic_binning</code> plotting app. Everything is drawn as the original draws it
(colors, color limits, edge widths, the zoom threshold, the cells, dots and labels); the differences are listed below.</p>
<table>
<tr><th>Element</th><th>Meaning</th></tr>
<tr><td>bin face color</td><td>mean anomaly (V &minus; V<sub>clim</sub>) of the profiles in the bin, at this depth level and in the
  time window; PRGn, cut so white is 0; limits = 5 % / 95 % quantiles of all values of the variable (all depth levels)</td></tr>
<tr><td>bin face opacity</td><td>sqrt(count / largest count), between 0.4 and 1</td></tr>
<tr><td>bin edge color</td><td>standard deviation (ddof 0) of the anomalies in the bin; cividis, 25 % to 95 % quantiles of the bins' std</td></tr>
<tr><td>bin edge width</td><td>0.01 pt per profile (none for one profile), capped so edges never cover the smallest bin</td></tr>
<tr><td>zoom level</td><td>sqrt(global area / view area); above 10 each bin with at most <i>max cells/bin</i> profiles splits into
  one cell per profile (Voronoi, 2 Lloyd steps), values sorted low to high from south to north; red dots are the true
  profile locations; numbers are the bin means</td></tr>
<tr><td>time</td><td>new: any window of whole years; "all years" is the whole record of the class</td></tr>
</table>
<p><b>Differences from the original.</b> (1) Only values the processing chain kept (weight &gt; 0) are shown; the original
used every non-missing value. (2) The color limits are computed per class over the whole record, so they stay fixed
while the time window changes. (3) The original showed one yearly file; this page shows any range of years, of all
files of the class. (4) If more than CELL_BUDGET cells would be in view, bins stay whole until you zoom further
(the browser's limit, not the original's). (5) Plain lon/lat map (as the original), drawn with WebGL (deck.gl).</p>
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
    html = html.replace('__PRGN__', lut('PRGn')).replace('__CIVIDIS__', lut('cividis'))
    html = html.replace('__GUIDE__', GUIDE_HTML.replace('CELL_BUDGET', '200,000'))
    (out / 'index.html').write_text(html)
    shutil.copytree(here / 'vendor', out / 'vendor', dirs_exist_ok=True)


def run(run_dir, out, sources=None, rebuild_bins=False):
    run_dir = common.yearly_dir(run_dir)
    out = Path(out)
    out.mkdir(parents=True, exist_ok=True)
    if rebuild_bins or not (out / 'bins.bin').exists():
        print('Bins (once):')
        build_bins(out)
    if not (out / 'land.json').exists():
        write_land(out)
    kept_bins = _kept_bins(out)
    groups = {}
    for f in sorted(run_dir.rglob('*.nc')):
        src, _ = common.classify(f, run_dir)
        if not sources or src in sources:
            groups.setdefault(src, []).append(f)
    tree, _ = _geodesic_tree()
    depth = None
    with netCDF4.Dataset(next(iter(groups.values()))[0]) as ds:
        depth = [float(x) for x in ds['prof_depth'][:]]
    index_path = out / 'index.json'
    index = json.loads(index_path.read_text()) if index_path.exists() else {'classes': {}}
    index.update({'run': str(run_dir), 'depth': depth, 'vars': list(VARS),
                  'units': {'T': '°C', 'S': 'psu'},
                  'created': time.strftime('%Y-%m-%dT%H:%M:%S')})
    print('Classes:')
    for src, files in sorted(groups.items()):
        index['classes'][src] = export_class(src, files, out, tree, kept_bins)
    index_path.write_text(json.dumps(index, indent=1))
    total = sum(p.stat().st_size for p in out.rglob('*') if p.is_file()) / 1e6
    write_page(out)
    print(f'Done: {out} ({total:.0f} MB)')
    return out


def upload(out, bucket, prefix, endpoint, profile='r2', dry_run=False):
    """Copy the folder to an S3-compatible bucket (Cloudflare R2) with the aws CLI.

    The data files are gzip bytes the page decompresses itself, so they are uploaded without Content-Encoding;
    the page is uploaded last, so it never points at data that isn't there yet.
    """
    import subprocess
    dest = f's3://{bucket}/{prefix.strip("/")}' if prefix else f's3://{bucket}'
    base = ['aws', 's3', 'sync', str(out), dest, '--endpoint-url', endpoint, '--profile', profile, '--only-show-errors']
    if dry_run:
        base.append('--dryrun')
    subprocess.run(base + ['--exclude', 'index.html', '--cache-control', 'public, max-age=86400'], check=True)
    subprocess.run(base + ['--exclude', '*', '--include', 'index.html', '--content-type', 'text/html; charset=utf-8',
                           '--cache-control', 'no-cache'], check=True)
    print(f'Uploaded {out} -> {dest}')


def cli(argv=None):
    p = argparse.ArgumentParser(description='Write the geodesic bin viewer data.')
    p.add_argument('--run', help='processed run folder (or its yearly/)')
    p.add_argument('--out', required=True, help='output folder (uploaded to the host as is)')
    p.add_argument('--sources', nargs='*', default=None, help='obs classes (default: all)')
    p.add_argument('--rebuild-bins', action='store_true', help='rebuild bins.bin even if it exists')
    p.add_argument('--page-only', action='store_true', help='only rewrite index.html (after editing viewer.html)')
    p.add_argument('--upload', metavar='BUCKET', help='then upload the folder to this R2 bucket')
    p.add_argument('--prefix', default='', help='folder in the bucket (e.g. run_20261027_142247)')
    p.add_argument('--endpoint', help='R2 endpoint, https://<account id>.r2.cloudflarestorage.com')
    p.add_argument('--profile', default='r2', help='aws CLI profile with the R2 API token (default: r2)')
    p.add_argument('--dry-run', action='store_true', help='show what the upload would copy')
    a = p.parse_args(argv)
    if a.page_only or (a.upload and not a.run):
        write_page(Path(a.out))
    else:
        run(a.run, a.out, a.sources, a.rebuild_bins)
    if a.upload:
        if not a.endpoint:
            p.error('--upload needs --endpoint')
        upload(Path(a.out), a.upload, a.prefix, a.endpoint, a.profile, a.dry_run)


if __name__ == '__main__':
    cli()
