"""
Data for the report's interactive map viewer (``profile-report --site``).

The viewer draws per-geodesic-bin maps in the browser, so any depth selection (single levels, 100 m
bands, groups of levels, flagged cost spikes), any 5-year time block and any cost part can be shown
without one image per combination. It needs two kinds of files, written under ``<site>/viewer/``:

- ``grid.json`` + ``grid_<view>.bin``: for each map view (global equal area, north and south polar
  stereographic) a raster of the geodesic bin id under each pixel (int16, -1 = land or off the map), and
  the coastline as polygons in the same pixel coordinates.
- ``<source>__<cost type>.bin`` (gzip): one record per non-empty (bin, level, time block) and variable,
  with the sum of the cost per datum for unsplit and anomaly, the point count, and per (bin, level) the
  split's mean term. ``index.json`` lists the files, depth levels, time blocks, groupings and flagged
  spike levels.

File layout (gzip of little-endian columns), per variable T then S: ``uint32 n``, then n values of each
column in turn: ``uint16 bin`` (1-based), ``uint8 level``, ``uint8 block``, ``float16`` mean cost per datum
unsplit, ``float16`` mean cost per datum of the split's anomaly term, ``uint16`` point count; then
``uint32 m`` and m values each of ``uint16 bin``, ``uint8 level``, ``float16`` mean term. float16 keeps
3 significant digits (values below ~6e-5 lose precision; they are ~0 on every color scale).
"""
import gzip
import json
from pathlib import Path

import netCDF4
import numpy as np

from ecco_insitu.diagnostics.profile_analyses import compute as pa

GROUPINGS = {'15': 15, '30': 30}       # the 97 levels grouped into about 15 or 30 consecutive groups


def _write_grid(out):
    """Pixel -> bin rasters and coastlines for the three views."""
    import cartopy.feature as cfeature
    import cartopy.crs as ccrs
    from ecco_insitu.diagnostics import common
    bl = common.BinLookup()
    views = {'global': (bl.projected_raster('global', npix=1400, proj=pa._global_proj()), pa._global_proj()),
             'NH': (bl.projected_raster('NH', npix=700, fill_box=pa._POLAR_EDGE_LAT), ccrs.NorthPolarStereo()),
             'SH': (bl.projected_raster('SH', npix=700, fill_box=pa._POLAR_EDGE_LAT), ccrs.SouthPolarStereo())}
    meta = {}
    land_geoms = list(cfeature.NaturalEarthFeature('physical', 'land', '110m').geometries())
    for name, ((raster, land, outside, extent), proj) in views.items():
        g = raster['a'].astype(np.int16)
        g[land | outside | (g == 0)] = -1
        ny, nx = g.shape
        (out / f'grid_{name}.bin').write_bytes(gzip.compress(g[::-1].tobytes(), 9))   # row 0 = top
        x0, x1, y0, y1 = extent
        polys = []
        for geom in land_geoms:
            pg = proj.project_geometry(geom, ccrs.PlateCarree())
            for part in getattr(pg, 'geoms', [pg]):
                if part.is_empty:
                    continue
                xy = np.asarray(part.exterior.coords)
                px = (xy[:, 0] - x0) / (x1 - x0) * nx
                py = (y1 - xy[:, 1]) / (y1 - y0) * ny
                keep = np.isfinite(px) & np.isfinite(py)
                if keep.sum() > 2:
                    polys.append(np.round(np.c_[px[keep], py[keep]], 1).ravel().tolist())
        meta[name] = {'nx': nx, 'ny': ny, 'coast': polys}
    (out / 'grid.json').write_text(json.dumps(meta))


def _groups(depth):
    """{'100m': [...], '15': [...], '30': [...]}: groups of level indices with labels. 100 m bands: levels
    with depth in [lo, lo + 100); bands with no level are skipped (below 2000 m the levels are 200 m apart).
    15 / 30 groups: the 97 levels split into that many consecutive groups."""
    d = np.asarray(depth, float)
    out = {'100m': []}
    for lo in np.arange(0, d.max() + 100, 100):
        lv = np.nonzero((d >= lo) & (d < lo + 100))[0].tolist()
        if lv:
            out['100m'].append({'label': f'{int(lo)}–{int(lo + 100)} m', 'levels': lv})
    for key, n in GROUPINGS.items():
        parts = np.array_split(np.arange(len(d)), n)
        out[key] = [{'label': f'{d[p[0]]:g}–{d[p[-1]]:g} m', 'levels': p.tolist()} for p in parts]
    return out


def _spikes(res, V, core, depth):
    """Depth levels where the mean cost per datum (true weights; unsplit or anomaly) stands out from the
    water column around it: at least 1.4 times the median of the levels 300-1000 m above and below it (at
    least 3 levels), and above 0.5. Broad bumps (several levels) are flagged at every level in them, so
    each flagged level is listed with its value and the reference median."""
    d = np.asarray(depth, float)
    out = []
    for part, sp in (('unsplit', 'unsplit'), ('anomaly', core)):
        r = res[('a02', V, sp, 'true')]
        m = np.where(r['count'] >= 30, r['mean'], np.nan)
        for k in range(len(m)):
            if not np.isfinite(m[k]) or m[k] <= 0.5:
                continue
            gap = np.abs(d - d[k])
            ref = m[(gap >= 300) & (gap <= 1000)]
            ref = ref[np.isfinite(ref)]
            if len(ref) >= 3 and m[k] >= 1.4 * np.median(ref):
                out.append({'level': k, 'part': part, 'value': float(m[k]), 'median': float(np.median(ref))})
    return out


def write(analysis_dir, out):
    """Write the viewer files for one profile-analyze folder into ``out`` (the site's viewer/ folder)."""
    analysis_dir, out = Path(analysis_dir), Path(out)
    out.mkdir(parents=True, exist_ok=True)
    meta = json.loads((analysis_dir / 'meta.json').read_text())
    core = 'split_own' if meta.get('core_split', pa.CORE_SPLIT) == 'own' else 'split'
    if not (out / 'grid.json').exists():
        _write_grid(out)
    index = {'files': {}, 'spikes': {}, 'core_split': meta.get('core_split', pa.CORE_SPLIT)}
    for src in meta['sources']:
        for ct in meta['cost_types']:
            f = sorted((analysis_dir / src).glob(f'analysis_{src}_*_{ct}.nc'))
            if not f:
                continue
            res = pa.load_results(f[-1])
            with netCDF4.Dataset(f[-1]) as nc:
                nc.set_auto_mask(False)
                if 'viewer__T__unsplit__sum' not in nc.variables:
                    continue
                if 'depth' not in index:
                    index['depth'] = [float(x) for x in nc['depth_m'][:]]
                    y0, nb = int(nc.viewer_year0), int(nc.viewer_block_years)
                    index['blocks'] = [f'{y0 + i * nb}–{y0 + i * nb + nb - 1}'
                                       for i in range(nc.dimensions['time_block'].size)]
                    index['groups'] = _groups(index['depth'])
                chunks = []
                for V in pa.VARS:
                    su = nc[f'viewer__{V}__unsplit__sum'][:]
                    sa = nc[f'viewer__{V}__{core}__sum'][:]
                    n = nc[f'viewer__{V}__unsplit__n'][:]
                    b, i, k = np.nonzero(n > 0)
                    nn = n[b, i, k]
                    cols = [(i + 1).astype('<u2'), k.astype('u1'), b.astype('u1'),
                            (su[b, i, k] / nn).astype('<f2'), (sa[b, i, k] / nn).astype('<f2'),
                            np.minimum(nn, 65535).astype('<u2')]
                    chunks += [np.uint32(len(b)).tobytes()] + [c.tobytes() for c in cols]
                    mt = nc[f'meanterm__{V}__true__by_bin_depth'][:] if index['core_split'] == 'own' and \
                        f'meanterm__{V}__true__by_bin_depth' in nc.variables else np.full((pa.N_BINS, 1), np.nan)
                    mi, mk = np.nonzero(np.isfinite(mt))
                    chunks += [np.uint32(len(mi)).tobytes(), (mi + 1).astype('<u2').tobytes(),
                               mk.astype('u1').tobytes(), mt[mi, mk].astype('<f2').tobytes()]
            name = f'{src}__{ct}.bin'
            (out / name).write_bytes(gzip.compress(b''.join(chunks), 6))
            index['files'][f'{src}|{ct}'] = name
            index['spikes'][f'{src}|{ct}'] = {V: _spikes(res, V, core, index['depth']) for V in pa.VARS}
    (out / 'index.json').write_text(json.dumps(index))
    return out
