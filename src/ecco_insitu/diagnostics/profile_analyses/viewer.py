"""
Data for the report's interactive map viewer (``profile-report --site``).

The viewer draws per-geodesic-bin maps in the browser, so any depth selection (single levels, 100 m
bands, groups of levels, flagged cost spikes), any 5-year time block and any cost part can be shown
without one image per combination. It needs two kinds of files, written under ``<site>/viewer/``:

- ``grid.json`` + ``grid_<view>.bin``: for each map view (global equal area, north and south polar
  stereographic) a raster of the geodesic bin id under each pixel (int16, -1 = land or off the map), and
  the coastline as polygons in the same pixel coordinates.
- ``<source>__<cost type>.bin`` (gzip, small, loaded first): per variable, one record per non-empty
  (bin, level) over all years, then the split's mean term per (bin, level).
- ``<source>__<cost type>__y<first>-<last>.bin`` (gzip, loaded when a time range needs them): the
  year-by-year data in pieces of YEARS_PER_PIECE years (year indices first..last): per variable, one record
  per non-empty (year, bin, level). Pieces keep every file under GitHub's 100 MB limit. ``index.json`` lists the files, depth levels, time blocks, groupings and flagged
  spike levels.

Record blocks (little-endian columns): ``uint32 n``, then n values of each column in turn: ``uint16 bin``
(1-based), ``uint8 level``, ``uint8 block`` (year index; 255 = all years), ``float16`` mean cost per datum
unsplit, ``float16`` mean cost per datum of the split's anomaly term, ``uint16`` point count, ``float16``
mean signed misfit V - V_ref (the bias diagnostic; degC or psu). The summary file has, per variable T then
S, one record block and then the mean term: ``uint32 m`` and m values each of ``uint16 bin``, ``uint8
level``, ``float16`` mean term. The years file has one record block per variable. float16 keeps
3 significant digits (values below ~6e-5 lose precision; they are ~0 on every color scale).
"""
import gzip
import json
from pathlib import Path

import netCDF4
import numpy as np

from ecco_insitu.diagnostics.profile_analyses import compute as pa

GROUPINGS = {'15': 15, '30': 30}       # the 97 levels grouped into about 15 or 30 consecutive groups


GUIDE_HTML = """
<p>The map explorer shows, for one obs class at a time, how well the observations fit a reference (climatology or
the model) in each <strong>geodesic bin</strong>: 10,242 cells of roughly equal area covering the globe. You choose
what to show with the rows of buttons; the maps redraw at once. Everything here is computed from the same saved
results as the report&rsquo;s figures.</p>

<h4>The rows of buttons</h4>
<table>
<tr><th>Row</th><th>Option</th><th>What it means</th></tr>
<tr><td rowspan="1">Source</td><td>CTD_WOD, PFL, &hellip;</td><td>The obs class (instrument type): CTD and bottle casts,
  gliders (GLD), ice-tethered profilers (ITP), seal tags (MEOP), moorings (MRB), Argo floats (PFL), XBTs. Each is
  analyzed on its own.</td></tr>
<tr><td rowspan="2">Cost type</td><td>vs_clim</td><td>Reference = the WOA23 climatology (1991&ndash;2020 average) at the profile&rsquo;s
  location, month and depths.</td></tr>
<tr><td>vs_model</td><td>Reference = ECCO v4r6, from its <em>monthly-mean</em> T and S sampled at the profile&rsquo;s grid cell,
  month and depths (approximate: not the model&rsquo;s value at the profile&rsquo;s exact time). No model values after 2025.</td></tr>
<tr><td>Variable</td><td>T, S</td><td>Temperature or salinity. S is grayed out for classes without salinity (XBT).</td></tr>
<tr><td rowspan="4">Cost</td><td>unsplit</td><td>The <strong>cost per datum</strong>: (&Delta;V/&sigma;)&sup2;, where &Delta;V is observation
  minus reference and &sigma; the assumed error of the observation. 1 means &ldquo;one assumed error away&rdquo;; 4 means two.
  Averaged over the bin&rsquo;s points in the selected depths and years. &ldquo;True weights&rdquo;: the area factor &gamma; and the spatial
  scaling factor sf are left out, so this is purely misfit relative to &sigma;.</td></tr>
<tr><td>anomaly (split)</td><td>The same, after removing the bin&rsquo;s <em>long-term mean misfit</em> at each depth level: it measures how well
  the <em>variations in time</em> are fit, ignoring a steady offset. This is the first term of ECCO&rsquo;s split cost
  (<code>cost_profiles.F</code>). The long-term mean is computed from this obs class&rsquo;s own points, all years.</td></tr>
<tr><td>mean term (split)</td><td>The second term of the split cost: the steady offset itself, (mean observation &minus; mean reference)&sup2;,
  weighted by the bin&rsquo;s combined error w&#772; = 1/max(s/n, s/9), s the mean of &sigma;&sup2; over the bin&rsquo;s n points (as
  <code>cost_profiles.F</code>). One value per bin and level; averaged over the selected levels. Defined over the whole record, so
  the time selection is off. Its values are typically much larger than the per-datum cost, since it doesn&rsquo;t shrink with
  the number of points: compare maps of the mean term with each other, not with the other parts.</td></tr>
<tr><td># points</td><td>The number of observations (profile values) in each bin for the selected depths and years: the
  &ldquo;n&rdquo; behind the averages. Bins with very few points give noisy averages.</td></tr>
<tr><td>Diagnostic</td><td>mean bias (obs &minus; reference)</td><td><strong>Not a cost term</strong>, so it has its own row. The
  average of the signed misfit, observation minus reference, over the bin&rsquo;s points in the selected depths and years, in
  &deg;C (T) or psu (S). Red: the observations are warmer or saltier than the reference; blue: colder or fresher. The color
  scale is symmetric about zero. Use it with the time range to see whether a steady offset is constant or
  <em>drifting</em> (model drift, a change in instruments or regions, or ocean change against a fixed climatology). It is
  related to the split&rsquo;s mean term (which is this offset, over the whole record, squared and weighted), but here it keeps
  its sign, has units, and can be limited to any years.</td></tr>
<tr><td rowspan="6">Depth</td><td>all depths</td><td>All depth levels together (each point counts once, so the upper ocean, which
  has more and closer-spaced levels, dominates).</td></tr>
<tr><td>single level</td><td>One of the 97 standard depth levels (2&nbsp;m to 6000&nbsp;m); pick it in the list.</td></tr>
<tr><td>100 m bands</td><td>All levels in a 100&nbsp;m range (0&ndash;100&nbsp;m, 100&ndash;200&nbsp;m, &hellip;). Below 2000&nbsp;m the levels are 200&nbsp;m
  apart, so some bands are skipped.</td></tr>
<tr><td>15 groups / 30 groups</td><td>The 97 levels split into 15 or 30 consecutive groups of about equal numbers of levels.</td></tr>
<tr><td>flagged spikes</td><td>Levels where the mean cost per datum of this class stands out: above 0.5 and at least 1.4 times the
  median of the levels 300&ndash;1000&nbsp;m above and below it (levels with at least 30 points). The list gives the level&rsquo;s
  value and that median, and whether it was flagged in the unsplit or anomaly cost. Use it to find <em>where</em> a peak in
  the report&rsquo;s cost-vs-depth plots (analyses 2 and 9) comes from.</td></tr>
<tr><td>(list)</td><td>The level, band, group or spike to show.</td></tr>
<tr><td rowspan="3">Time</td><td>all years</td><td>Every year of the record. Fast: uses a small all-years summary.</td></tr>
<tr><td>the two handles</td><td>A range of years: drag the left handle for the first year and the right one for the last. The map
  shows the average over every point in those years. Choosing years downloads that class&rsquo;s year-by-year data in 5-year
  pieces, only those the range needs (up to about 20&nbsp;MB per piece for Argo, a few MB for most classes), so the first
  range in a period can take a few seconds; after that, changes are quick.</td></tr>
<tr><td>single year</td><td>Locks the handles together, so dragging moves a one-year window.</td></tr>
</table>

<h4>The maps and the color bar</h4>
<p><strong>Maps:</strong> global (equal-area projection, starting at the western tip of Australia), and north and south polar
stereographic views. White = no data in that bin for the selection; gray = land. Hover a bin to see its number and value.</p>
<p><strong>Color bar:</strong> ticks give the values; by default the scale runs from 0 to the 98th percentile of the bins shown, and
values above it get the top color (so a few extreme bins don&rsquo;t wash out the rest). <strong>max</strong>: set the top of the scale
yourself, e.g. to compare two selections on the same scale (the scale otherwise changes with each selection).
<strong>auto</strong>: back to the 98th percentile. <strong>log color scale</strong>: spreads small values out (covers three powers of
ten below the max).</p>

<h4>Reading it</h4>
<ul>
<li>Values near 1 mean the misfit is about what the assumed errors &sigma; allow; much larger values point to a poor fit, a
too-small &sigma;, or bad data. With few points (# points) a single bad profile can dominate a bin.</li>
<li>If unsplit is high but anomaly is low in a region, the misfit there is mostly a steady offset (see the mean term).
  The mean bias shows that offset&rsquo;s sign and size; compare early and late year ranges to see if it changes.</li>
<li>The time range shows whether a high cost is persistent or limited to some years (e.g. one cruise or float).</li>
<li>Data notes, open decisions and corrections are on the report&rsquo;s <em>Notes &amp; errata</em> page.</li>
</ul>
<p style="color:#777">Precision: values are stored to 3 significant digits. Time resolution: one year. Bins are the same geodesic bins
as the report and the ECCO split cost.</p>
"""


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


YEARS_PER_PIECE = 5      # year-by-year data is written in pieces of this many years

_COLS = (('bin', '<u2'), ('lev', 'u1'), ('blk', 'u1'), ('mu', '<f2'), ('ma', '<f2'), ('n', '<u2'), ('bias', '<f2'))


def _cells(nc, V, core, nd):
    """Per non-empty (year, bin, level): mean cost per datum unsplit and anomaly, point count, mean bias."""
    if f'viewer__{V}__unsplit__cell' not in nc.variables:
        return None
    cell = nc[f'viewer__{V}__unsplit__cell'][:]
    su = nc[f'viewer__{V}__unsplit__sum'][:].astype(np.float64)
    nn = nc[f'viewer__{V}__unsplit__n'][:].astype(np.float64)
    ca = nc[f'viewer__{V}__{core}__cell'][:]          # a subset of unsplit's cells (points with a bin mean)
    sa = np.full(len(cell), np.nan); na = np.zeros(len(cell))
    pos = np.searchsorted(cell, ca)
    sa[pos], na[pos] = nc[f'viewer__{V}__{core}__sum'][:], nc[f'viewer__{V}__{core}__n'][:]
    mis = nc[f'viewer__{V}__unsplit__misfit_sum'][:].astype(np.float64) \
        if f'viewer__{V}__unsplit__misfit_sum' in nc.variables else np.full(len(cell), np.nan)
    with np.errstate(invalid='ignore', divide='ignore'):
        return {'bin': (cell // nd) % (pa.N_BINS + 1), 'lev': cell % nd, 'blk': cell // nd // (pa.N_BINS + 1),
                'mu': su / nn, 'ma': np.where(na > 0, sa / na, np.nan), 'n': nn, 'na': na, 'bias': mis / nn}


def _all_years(r):
    """The same, summed over the years (block = 255 marks 'all years')."""
    if not len(r['bin']):
        return r
    key = r['bin'].astype(np.int64) * 256 + r['lev']
    u, inv = np.unique(key, return_inverse=True)
    n = np.bincount(inv, weights=r['n'])
    na = np.bincount(inv, weights=r['na'])
    ok = np.isfinite(r['ma'])
    with np.errstate(invalid='ignore', divide='ignore'):
        return {'bin': u // 256, 'lev': u % 256, 'blk': np.full(len(u), 255),
                'mu': np.bincount(inv, weights=r['mu'] * r['n']) / n,
                'ma': np.where(na > 0, np.bincount(inv, weights=np.where(ok, r['ma'] * r['na'], 0)) / na, np.nan),
                'n': n, 'bias': np.bincount(inv, weights=np.nan_to_num(r['bias']) * r['n']) / n}


def _pack(r):
    """uint32 count, then each column of _COLS."""
    cols = [np.asarray(r[c]) for c, _ in _COLS]
    cols[5] = np.minimum(cols[5], 65535)
    return [np.uint32(len(cols[0])).tobytes()] + [c.astype(dt).tobytes() for c, (_, dt) in zip(cols, _COLS)]


def _pack_meanterm(nc, V, core_split):
    name = f'meanterm__{V}__true__by_bin_depth'
    mt = nc[name][:] if core_split == 'own' and name in nc.variables else np.full((pa.N_BINS, 1), np.nan)
    mi, mk = np.nonzero(np.isfinite(mt))
    return [np.uint32(len(mi)).tobytes(), (mi + 1).astype('<u2').tobytes(), mk.astype('u1').tobytes(),
            mt[mi, mk].astype('<f2').tobytes()]


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
                if 'viewer__T__unsplit__cell' not in nc.variables:
                    continue
                nd = nc.dimensions['depth'].size
                if 'depth' not in index:
                    index['depth'] = [float(x) for x in nc['depth_m'][:]]
                    y0, nb, nblk = int(nc.viewer_year0), int(nc.viewer_block_years), int(nc.viewer_n_blocks)
                    index['blocks'] = [str(y0 + i * nb) if nb == 1 else f'{y0 + i * nb}–{y0 + i * nb + nb - 1}'
                                       for i in range(nblk)]
                    index['groups'] = _groups(index['depth'])
                summary, recs = [], {}
                for V in pa.VARS:
                    rec = _cells(nc, V, core, nd)
                    if rec is None:                      # no data for this variable (e.g. S for XBT)
                        rec = {c: np.zeros(0, dt) for c, dt in _COLS}
                    recs[V] = rec
                    # all years: one record per (bin, level); 'n'-weighted means
                    summary += _pack(_all_years(rec))
                    summary += _pack_meanterm(nc, V, index['core_split'])
            name = f'{src}__{ct}'
            (out / f'{name}.bin').write_bytes(gzip.compress(b''.join(summary), 6))
            index['files'][f'{src}|{ct}'] = f'{name}.bin'
            # year-by-year data in pieces of YEARS_PER_PIECE years (GitHub allows at most 100 MB per file,
            # and the browser then downloads only the pieces a time range needs)
            pieces = []
            nblk = len(index['blocks'])
            for p0 in range(0, nblk, YEARS_PER_PIECE):
                p1 = min(p0 + YEARS_PER_PIECE, nblk) - 1
                chunk = []
                for V in pa.VARS:
                    r = recs[V]
                    m = (r['blk'] >= p0) & (r['blk'] <= p1) if len(r['blk']) else np.zeros(0, bool)
                    chunk += _pack({c: np.asarray(r[c])[m] for c, _ in _COLS})
                fname = f'{name}__y{p0:02d}-{p1:02d}.bin'
                (out / fname).write_bytes(gzip.compress(b''.join(chunk), 6))
                pieces.append({'file': fname, 'first': p0, 'last': p1})
            index.setdefault('years_files', {})[f'{src}|{ct}'] = pieces
            index['spikes'][f'{src}|{ct}'] = {V: _spikes(res, V, core, index['depth']) for V in pa.VARS}
    (out / 'index.json').write_text(json.dumps(index))
    return out
