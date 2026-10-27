"""
**Work in progress** (see :mod:`ecco_insitu.model_comparison`). Reproduces the v4r6 iter47
per-type cost bars; not yet checked against a model run's own cost output.

Offline profile cost, computed as ECCO v4r6's cost_profiles.F does it (ECCO-GROUP/ECCO-v4-Configurations,
"ECCOv4 Release 6/code/cost_profiles.F", with ALLOW_PROFILES_SAMPLESPLIT_COST and
ALLOW_PROFILES_CLIMMASK defined, as in that release's PROFILES_OPTIONS.h).

Per variable V (T, S, ...), with  d = prof_V,  m = model equivalent (prof_Vestim),
w = prof_Vweight, and the averaging bin b = prof_bin_id_a (truncated to an integer)::

  CLIMMASK      w = 0 where prof_Vclim <= -990                       (profiles_readvector.F)
  model mask    mask = 1 where m != 0, else 0                         (no model value is 0)
  bin means     over all files, per (bin, depth value), using points with w > 0 and mask = 1:
                  m_mean = mean(m),  d_mean = mean(d),  s = mean(1/w),  n = count
                  w_mean = 1 / max(s / n, s / profiles_mean_indsamples)        (default 9)
                Depths are matched across files by exact value (prof_depth_comb in the model).
  anomaly cost  sum over points with w > 0 and d_mean != 0 of
                  w * mask * (m - d - m_mean + d_mean)^2,     count = sum of mask
                (reported per file and variable, as 'cost_profiles(file,var)')
  mean cost     sum over (bin, depth) with w_mean > 0, d_mean != 0 and m_mean != 0 of
                  w_mean * (m_mean - d_mean)^2,               count = number of terms
                (reported per variable, as 'cost_profiles_mean(var)')
  total         sum of mult_profiles * anomaly cost + mult_profiles_mean * mean cost
                (all multipliers are 1 in v4r6's data.profiles)

**Versions with weight factors removed** (to see what drives a difference in cost). The weight
is w = gamma * sf_V / sigma^2, where gamma = ``prof_area_gamma`` and sf_V is the spatial scaling
("decimation") factor: ``prof_<V>spatial_scaling_factor`` in this package's files,
``prof_spatial_scaling_factor`` in older ones such as ECCO v4r6's. Besides ``final`` (the model's
cost), the anomaly cost is also reported as

- ``no_sf``: each point's weight divided by sf_V (spatial scaling removed);
- ``no_gamma_sf``: divided by gamma * sf_V, i.e. (misfit / sigma)^2.

The bin means are the model's in all three versions; only each point's weight changes.

**Breakdowns:** by variable, file (so dataset type), year, latitude band (``LAT_EDGES``) and depth band
(``DEPTH_EDGES``), in ``cost_breakdown.csv`` (one row per group, with the point count and the cost of
each version). Group and plot it any way you like, e.g. with ``ecco-insitu model-cost-bars``.

Model values are read from the data file if it has prof_Vestim, else from <name>_model.nc in
--estim (e.g. made by ``ecco-insitu mds-to-nc``, or MITprof_gcm2nc.m).
-9999 in prof_Vestim means "no model value" and is treated as 0, as in the model's binary output.

Usage::

    ecco-insitu model-cost --data <dir with profile .nc> [--estim <dir with *_model.nc>]
                           [--vars T S] [--indsamples 9] [--out <dir>] [-n 8]

"""
import argparse
import sys
import glob
import json
import os
from concurrent.futures import ProcessPoolExecutor

import netCDF4
import numpy as np

CHUNK = 20000            # profiles read at a time
NAVGBIN = 10246          # bin ids 1..10245 (NAVGBINMAX in v4r6's PROFILES_SIZE.h); 0 unused


def estim_path(path, estim_dir):
    """File holding the model values for data file ``path``: ``path`` itself if ``estim_dir`` is
    None, else ``<estim_dir>/<name>_model.nc``."""
    if estim_dir is None:
        return path
    name = os.path.basename(path)[:-3]
    return os.path.join(estim_dir, name + '_model.nc')


def _factor(ds, names, i0, i1):
    """The first of ``names`` present in ``ds``, for profiles i0:i1, as floats; None if none is present.
    Missing values (fill, non-finite, <= 0) become NaN."""
    for n in names:
        if n in ds.variables:
            x = np.asarray(ds[n][i0:i1], float)
            return np.where(np.isfinite(x) & (x > 0) & (x < 1e30), x, np.nan)
    return None


def read_chunk(path, epath, V, i0, i1, extras=False):
    """d, m, w (after CLIMMASK) for profiles i0:i1, and bin ids (int, -1 where unusable).

    With ``extras``, also a dict of per-profile values for the breakdowns: ``gamma``
    (``prof_area_gamma``), ``sf`` (the variable's spatial scaling factor:
    ``prof_<V>spatial_scaling_factor`` in this package's files, ``prof_spatial_scaling_factor`` in
    older ones), ``lat`` and ``year``. A factor missing from the file is NaN, and the versions that
    need it are then not computed for that file.
    """
    with netCDF4.Dataset(path) as ds:
        ds.set_auto_mask(False)
        d = np.asarray(ds[f'prof_{V}'][i0:i1], float)
        w = np.asarray(ds[f'prof_{V}weight'][i0:i1], float)
        if f'prof_{V}clim' in ds.variables:
            c = np.asarray(ds[f'prof_{V}clim'][i0:i1], float)
            w = np.where(c <= -990, 0.0, w)
        b = np.asarray(ds['prof_bin_id_a'][i0:i1], float)
        if extras:
            n = i1 - i0
            nan = np.full(n, np.nan)
            g = _factor(ds, ['prof_area_gamma'], i0, i1)
            sf = _factor(ds, [f'prof_{V}spatial_scaling_factor', 'prof_spatial_scaling_factor'], i0, i1)
            ex = {'gamma': nan if g is None else g, 'sf': nan if sf is None else sf,
                  'lat': np.asarray(ds['prof_lat'][i0:i1], float),
                  'year': (np.asarray(ds['prof_YYYYMMDD'][i0:i1], float) // 10000).astype(np.int64)}
    with netCDF4.Dataset(epath) as ds:
        ds.set_auto_mask(False)
        m = np.asarray(ds[f'prof_{V}estim'][i0:i1], float)
    m = np.where(m == -9999, 0.0, m)
    with np.errstate(invalid='ignore'):
        bi = np.where(np.isfinite(b) & (b >= 1) & (b < NAVGBIN), np.trunc(np.nan_to_num(b)), -1).astype(np.int64)
    if extras:
        return d, m, w, bi, ex
    return d, m, w, bi


def pass1(args):
    """Bin sums for one chunk: sum m, sum d, sum 1/w, count, per flat (bin, comb level)."""
    path, epath, V, i0, i1, lev, nlev = args
    d, m, w, bi = read_chunk(path, epath, V, i0, i1)
    with np.errstate(invalid='ignore'):
        use = (w > 0) & (m != 0) & (bi[:, None] >= 1)
    fin = np.isfinite(d) & np.isfinite(m) & np.isfinite(w)
    ok = use & fin
    flat = (bi[:, None] * nlev + lev[None, :])[ok]
    size = NAVGBIN * nlev
    bad = int((use & ~fin).sum())
    return (np.bincount(flat, weights=m[ok], minlength=size),
            np.bincount(flat, weights=d[ok], minlength=size),
            np.bincount(flat, weights=1.0 / w[ok], minlength=size),
            np.bincount(flat, minlength=size).astype(float), bad)


#: The three versions of the cost, by which weight factors are removed (see the module description).
VERSIONS = ('final', 'no_sf', 'no_gamma_sf')
#: Latitude bands (edges, degrees) and depth bands (edges, m) of the breakdown.
LAT_EDGES = (-90, -60, -30, 0, 30, 60, 90)
DEPTH_EDGES = (0, 200, 1000, 2000, 7000)


def _band(x, edges):
    """Index of the band each value falls in (``edges[i] <= x < edges[i+1]``; the top edge included)."""
    return np.clip(np.searchsorted(edges, x, side='right') - 1, 0, len(edges) - 2)


def pass2(args):
    """Anomaly cost for one chunk, given the bin means: totals, and a breakdown table.

    Returns ``(cost, count, table)``. ``cost`` and ``count`` are the model's numbers (version
    ``final``). ``table`` maps ``(year, lat band, depth band)`` to ``[count, final, no_sf,
    no_gamma_sf, n_with_factors]``. The two reduced versions divide each point's weight by its
    spatial factor, and by gamma times the spatial factor. They are summed only over points where
    those factors are known, counted in ``n_with_factors``.
    """
    path, epath, V, i0, i1, lev, nlev, mmean, dmean, depth = args
    d, m, w, bi, ex = read_chunk(path, epath, V, i0, i1, extras=True)
    flat = np.where(bi[:, None] >= 1, bi[:, None] * nlev + lev[None, :], 0)
    mm, dm = mmean[flat], dmean[flat]            # flat 0 is bin 0, whose means are 0
    mask = (m != 0).astype(float)
    with np.errstate(invalid='ignore'):
        sel = (w > 0) & (dm != 0) & np.isfinite(d) & np.isfinite(m) & np.isfinite(w)
    r2 = (m - d - mm + dm) ** 2
    final = np.where(sel, w * mask * r2, 0.0)
    cnt = np.where(sel, mask, 0.0)
    g = ex['gamma'][:, None]
    sf = ex['sf'][:, None]
    with np.errstate(invalid='ignore', divide='ignore'):
        no_sf = final / sf
        no_gs = final / (g * sf)
    has = sel & np.isfinite(no_gs) & (mask > 0)
    no_sf = np.where(has, no_sf, 0.0)
    no_gs = np.where(has, no_gs, 0.0)

    shape = final.shape
    yr = np.broadcast_to(ex['year'][:, None], shape)
    lb = np.broadcast_to(_band(ex['lat'], LAT_EDGES)[:, None], shape)
    db = np.broadcast_to(_band(np.asarray(depth, float), DEPTH_EDGES)[None, :], shape)
    table = {}
    pick = cnt > 0
    keys = np.stack([yr[pick], lb[pick], db[pick]], axis=1)
    vals = np.stack([cnt[pick], final[pick], no_sf[pick], no_gs[pick], has[pick].astype(float)], axis=1)
    if keys.size:
        uk, inv = np.unique(keys, axis=0, return_inverse=True)
        sums = np.zeros((len(uk), vals.shape[1]))
        np.add.at(sums, inv.ravel(), vals)
        table = {tuple(int(x) for x in k): sums[i].tolist() for i, k in enumerate(uk)}
    return float(final.sum()), float(cnt.sum()), table


#: Dataset type from the file path; checked in this order (MEOP files also contain 'PFL').
TYPE_RULES = (('MEOP', 'MEOP'), ('ITP', 'ITP'), ('LNW', 'LNW'), ('XBT', 'XBT'), ('GLD', 'GLD'),
              ('MRB', 'MRB'), ('CTD', 'CTD_OSD'), ('PFL', 'PFL'))
BREAKDOWN_COLUMNS = ('var', 'type', 'file', 'year', 'lat_band', 'depth_band', 'count',
                     'cost_final', 'cost_no_sf', 'cost_no_gamma_sf', 'count_with_factors')


def dataset_type(path):
    """Dataset type of a file path: PFL, MEOP, CTD_OSD, ITP, MRB, GLD, XBT, LNW, or 'other'."""
    up = path.upper()
    for key, name in TYPE_RULES:
        if key in up:
            return name
    return 'other'


def _band_label(i, edges, unit):
    return f'{edges[i]}..{edges[i + 1]}{unit}'


def breakdown_rows(breakdown):
    """Flatten the breakdown into rows (dicts with :data:`BREAKDOWN_COLUMNS`)."""
    rows = []
    for (V, f), tab in sorted(breakdown.items()):
        for (yr, lb, db), (n, cf, cs, cg, nf) in sorted(tab.items()):
            rows.append({'var': V, 'type': dataset_type(f), 'file': f, 'year': yr,
                         'lat_band': _band_label(lb, LAT_EDGES, ''), 'depth_band': _band_label(db, DEPTH_EDGES, 'm'),
                         'count': n, 'cost_final': cf, 'cost_no_sf': cs, 'cost_no_gamma_sf': cg,
                         'count_with_factors': nf})
    return rows


def write_breakdown_csv(rows, path):
    """Write the breakdown rows to a CSV file: one row per variable, file, year, latitude band and
    depth band, with the point count and the anomaly cost of each version. Cost per datum is
    ``cost_final / count``; for the reduced versions, divide by ``count_with_factors``."""
    import csv
    with open(path, 'w', newline='') as fh:
        wr = csv.DictWriter(fh, fieldnames=BREAKDOWN_COLUMNS)
        wr.writeheader()
        wr.writerows(rows)


def summarize(rows, by=('var', 'type')):
    """Cost per datum of each version, grouped by the columns ``by``.

    Returns
    -------
    dict
        ``{group key tuple: {'count': n, 'final': ..., 'no_sf': ..., 'no_gamma_sf': ...}}``, the
        per-datum values (the reduced versions over the points whose factors are known).
    """
    acc = {}
    for r in rows:
        k = tuple(r[c] for c in by)
        a = acc.setdefault(k, [0.0] * 5)
        for i, c in enumerate(('count', 'cost_final', 'cost_no_sf', 'cost_no_gamma_sf', 'count_with_factors')):
            a[i] += r[c]
    out = {}
    for k, (n, cf, cs, cg, nf) in acc.items():
        out[k] = {'count': n, 'final': cf / n if n else float('nan'),
                  'no_sf': cs / nf if nf else float('nan'), 'no_gamma_sf': cg / nf if nf else float('nan')}
    return out


def print_type_summary(rows):
    """Print the anomaly cost per datum by variable and dataset type, for the three versions."""
    s = summarize(rows)
    print('\nAnomaly cost per datum by type (final = what the model sees; no_sf = spatial scaling removed;'
          ' no_gamma_sf = gamma and spatial scaling removed):')
    print(f'  {"var":3s} {"type":8s} {"count":>12s} {"final":>9s} {"no_sf":>9s} {"no_gamma_sf":>12s}')
    for (V, t), v in sorted(s.items()):
        print(f'  {V:3s} {t:8s} {v["count"]:12.0f} {v["final"]:9.4f} {v["no_sf"]:9.4f} {v["no_gamma_sf"]:12.4f}')


def cli(argv=None):
    """Command-line entry point (``ecco-insitu`` subcommand)."""
    p = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    p.add_argument('--data', required=True, help='folder with the profile files (searched recursively)')
    p.add_argument('--estim', help='folder with <name>_model.nc files holding prof_Vestim (if the data files lack it)')
    p.add_argument('--vars', nargs='+', default=['T', 'S'])
    p.add_argument('--indsamples', type=float, default=9.0, help='profiles_mean_indsamples (v4r6: 9)')
    p.add_argument('--out', help='folder for cost_summary.json (default: print only)')
    p.add_argument('-n', type=int, default=4, help='parallel workers (default 4)')
    a = p.parse_args(argv)

    files = sorted(f for f in glob.glob(os.path.join(a.data, '**', '*.nc'), recursive=True)
                   if not f.endswith('_model.nc'))
    # per file: profile count, depth levels, and which variables it has
    info = {}
    for f in files:
        ep = estim_path(f, a.estim)
        if not os.path.exists(ep):
            print(f'{os.path.basename(f)}: no model file {ep}, skipped'); continue
        with netCDF4.Dataset(f) as ds, netCDF4.Dataset(ep) as es:
            vs = [V for V in a.vars if f'prof_{V}' in ds.variables and f'prof_{V}weight' in ds.variables
                  and f'prof_{V}estim' in es.variables]
            if 'prof_bin_id_a' not in ds.variables:
                print(f'{os.path.basename(f)}: no prof_bin_id_a, skipped'); continue
            info[f] = (ep, len(ds.dimensions['iPROF']), np.asarray(ds['prof_depth'][:], float), vs)
    comb = np.unique(np.concatenate([v[2] for v in info.values()]))      # combined depths, exact values
    comb = comb[comb >= 0]
    nlev = len(comb)
    lev_of = {f: np.searchsorted(comb, v[2]) for f, v in info.items()}

    result = {'settings': {'indsamples': a.indsamples, 'files': len(info), 'combined_levels': nlev}}
    breakdown = {}                     # (var, file) -> {(year, lat band, depth band): [count, final, no_sf, no_gs, n_f]}
    total = 0.0
    with ProcessPoolExecutor(a.n) as ex:
        for V in a.vars:
            fl = [f for f in info if V in info[f][3]]
            jobs = [(f, info[f][0], V, i0, min(i0 + CHUNK, info[f][1]), lev_of[f], nlev)
                    for f in fl for i0 in range(0, info[f][1], CHUNK)]
            if not jobs:
                continue
            Sm = Sd = Ss = N = 0; bad = 0
            for sm, sd, ss, n, bd in ex.map(pass1, jobs):
                Sm, Sd, Ss, N, bad = Sm + sm, Sd + sd, Ss + ss, N + n, bad + bd
            with np.errstate(invalid='ignore', divide='ignore'):
                mmean = np.where(N > 0, Sm / N, 0.0); dmean = np.where(N > 0, Sd / N, 0.0)
                s = np.where(N > 0, Ss / N, 0.0)
                sig2 = np.maximum(s / np.where(N > 0, N, 1), s / a.indsamples)
                wmean = np.where((N > 0) & (sig2 != 0), 1.0 / np.where(sig2 != 0, sig2, 1), 0.0)
            sel = (wmean > 0) & (dmean != 0) & (mmean != 0)
            jmean = float(np.sum(wmean[sel] * (mmean[sel] - dmean[sel]) ** 2)); nmean = int(sel.sum())

            per_file = {}
            jobs2 = [j + (mmean, dmean, info[j[0]][2]) for j in jobs]
            for j, (c, n, table) in zip(jobs2, ex.map(pass2, jobs2)):
                k = os.path.relpath(j[0], a.data)
                pc, pn = per_file.get(k, (0.0, 0.0)); per_file[k] = (pc + c, pn + n)
                ftab = breakdown.setdefault((V, k), {})
                for key, vals in table.items():
                    acc = ftab.setdefault(key, [0.0] * 5)
                    for i, x in enumerate(vals):
                        acc[i] += x
            janom = sum(c for c, _ in per_file.values()); nanom = sum(n for _, n in per_file.values())
            total += janom + jmean
            result[V] = {'anomaly_cost': janom, 'anomaly_count': nanom, 'mean_cost': jmean, 'mean_count': nmean,
                         'per_file': {k: {'cost': c, 'count': n} for k, (c, n) in per_file.items()},
                         'nonfinite_points_used': bad}

            print(f'\n{V}:')
            for i, (k, (c, n)) in enumerate(sorted(per_file.items()), 1):
                print(f'  cost_profiles({i:2d},{V})= {c:12.5E} {n:12.5E}   {k}')
            print(f'  cost_profiles_mean({V})= {jmean:12.5E} {float(nmean):12.5E}')
            print(f'  anomaly total {janom:.6E} (count {nanom:.0f}), mean {jmean:.6E}; cost per anomaly point '
                  f'{janom / nanom if nanom else float("nan"):.4f}')
            if bad:
                print(f'  WARNING: {bad} points with weight > 0 have non-finite data, model or weight values; left out')
    result['total'] = total
    print(f'\nTotal profiles cost (all multipliers 1): {total:.6E}')
    rows = breakdown_rows(breakdown)
    print_type_summary(rows)
    if a.out:
        os.makedirs(a.out, exist_ok=True)
        with open(os.path.join(a.out, 'cost_summary.json'), 'w') as fh:
            json.dump(result, fh, indent=1)
        write_breakdown_csv(rows, os.path.join(a.out, 'cost_breakdown.csv'))
        print('wrote', os.path.join(a.out, 'cost_summary.json'), 'and cost_breakdown.csv')


if __name__ == '__main__':
    sys.exit(cli())
