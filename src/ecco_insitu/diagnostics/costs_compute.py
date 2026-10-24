"""
Recompute model-data cost vs climatology from the FINAL chain files, store it
per profile with provenance, and collect statistics per source and year.

    cost_V = (prof_V - prof_Vclim)^2 * prof_Vweight        V in {T, S}

computed at every level where V, Vclim and Vweight are all real values (not
fill) and Vweight > 0, i.e. points that actually enter the fit.

The final weight is  gamma / sigma^2 * sf_V : sigma from the CTD sigma files
(step 04), gamma = prof_area_gamma (grid-cell area factor, step 10 since
2026-10-12 and step 05 before; shared by T
and S), sf_V = prof_Vspatial_scaling_factor (step 11, 1/n, computed separately
for T and S; with a single shared factor sf_T = sf_S and nothing here changes).
Three versions of the cost (see KINDS)::

  final  (a)  nothing removed: weight as written — what ECCO sees
              = gamma * sf_V * (dV/sigma)^2
  nosf   (b)  step 11 sf_V removed, gamma kept = gamma * (dV/sigma)^2
              (what QC tested in runs where gamma came before QC)
  true   (c)  gamma and sf_V removed = (dV/sigma)^2, the normalized misfit
              (what QC tests now that it runs before gamma)

Fill handling::

  T, S, Tclim, Sclim : non-finite or |x| >= 900 (no physical T/S is near that;
                       catches -9999, 9999, 99999, 9.97e36, 999, ...)
  weights            : non-finite, |x| >= 1e30, or |x| exactly 999, 9999,
                       99999, ... (real weights can exceed 9000 when sigma is
                       small, so they are kept; how many is reported)

Outputs (in <out_dir>)::

  costs/<source>/<file stem>__cost.nc   per file: cost_T, cost_S (iPROF, iDEPTH,
                                        final weights, NaN where excluded), the
                                        scaling factors, and provenance (descr,
                                        lat, lon, date, depth, source, file)
  profiles.parquet                      one row per profile: provenance, gamma,
                                        sfV, and per V: n used, depth of max,
                                        mean / max for each version
                                        (V_mean_final, V_max_true, ...)
  stats_source_year.csv                 per source x year x V: n, and min, max,
                                        mean, median for each version
  stats_source.csv                      same, all years
  fills.csv                             per source: points dropped as fill
  outliers.csv                          top points by true cost per V, with all
                                        three versions

Run::

    ecco-insitu report-costs            (computes, then writes the report)
    python -m ecco_insitu.diagnostics.costs_compute [--processed RUN] [--out DIR] [-n 4]   (compute only)
"""

import argparse
import concurrent.futures
import datetime
import json
import os
import sys
import time
from collections import defaultdict
from pathlib import Path

import netCDF4
import numpy as np
import pandas as pd

import ecco_insitu.diagnostics.common as common  # noqa: E402  (run selection, source/subtype parsing)

VARS = ('T', 'S')
REPORT_PARENT = Path(common._CFG['diagnostics']['cost_report_dir'] or '.')   # from the chain settings
N_OUTLIERS = 200
_FILL_9S = np.array([999, 9999, 99999, 999999, 9999999], dtype=np.float64)

# The three cost versions, in report order.  Each divides the final cost by
# what it removes: 'none', 'sf' (step 11) or 'gamma*sf' (step 05 and 11).
KINDS = ('final', 'nosf', 'true')
# Labels name the variable's own step 11 factor (sf_T, sf_S); V='V' for text
# that covers both.  fmt: 'text' sf_T, 'html' sf<sub>T</sub>, 'mpl' $sf_T$.
_KIND_LABEL = {'final': 'a) nothing removed (final weight)',
               'nosf':  'b) step 11 {sf} removed (γ kept)',
               'true':  'c) γ and {sf} removed = (ΔV/σ)²'}
_KIND_SHORT = {'final': 'a) final', 'nosf': 'b) no {sf}', 'true': 'c) true'}
_SF_FMT = {'text': 'sf_{V}', 'html': 'sf<sub>{V}</sub>', 'mpl': '$sf_{{{V}}}$'}


def sf_name(V='V', fmt='text'):
    """Name of the step 11 scaling factor for variable ``V`` as text, HTML or matplotlib math
    (``fmt`` = ``'text'``, ``'html'`` or ``'mpl'``).
    """
    return _SF_FMT[fmt].format(V=V)


def kind_label(k, V='V', fmt='text'):
    """Long label of cost version ``k`` (``'final'``, ``'nosf'``, ``'true'``) for variable ``V``."""
    return _KIND_LABEL[k].format(sf=sf_name(V, fmt))


def kind_short(k, V='V', fmt='text'):
    """Short label (``a) final``, ``b) no sf``, ``c) true``) of cost version ``k`` for variable ``V``."""
    return _KIND_SHORT[k].format(sf=sf_name(V, fmt))


def _read(v, name):
    if name not in v:
        return None
    x = v[name][:]
    return np.asarray(np.ma.filled(np.ma.asarray(x, dtype=np.float64), np.nan))


def _is_fill_value(x):
    """For T, S, clim."""
    return ~np.isfinite(x) | (np.abs(x) >= 900)


def _is_fill_weight(w):
    a = np.abs(w)
    return ~np.isfinite(w) | (a >= 1e30) | np.isin(a, _FILL_9S)


def _descr(v, n):
    if 'prof_descr' not in v:
        return np.full(n, '', dtype=object)
    d = netCDF4.chartostring(np.ma.filled(v['prof_descr'][:], b''))
    return np.char.strip(d.astype(str)).astype(object)


def _summ(x):
    if x.size == 0:
        return dict(n=0, min=np.nan, max=np.nan, mean=np.nan, median=np.nan)
    return dict(n=int(x.size), min=float(x.min()), max=float(x.max()),
                mean=float(x.mean()), median=float(np.median(x)))


def _process_file(path, run_dir, out_dir):
    """Return (per-profile DataFrame, {V: {kind: used-point array}},
    fill counts, outlier rows)."""
    source, subtype = common.classify(path, run_dir)
    with netCDF4.Dataset(path) as nc:
        nc.set_auto_mask(False)
        v = nc.variables
        n = nc.dimensions['iPROF'].size
        depth = _read(v, 'prof_depth')
        lat, lon = _read(v, 'prof_lat'), _read(v, 'prof_lon')
        ymd = _read(v, 'prof_YYYYMMDD').astype(np.int64)
        hms = _read(v, 'prof_HHMMSS').astype(np.int64)
        descr = _descr(v, n)
        gamma = _read(v, 'prof_area_gamma')
        bin_a = _read(v, 'prof_bin_id_a')

        prof = {'source': source, 'subtype': subtype, 'file': path.name,
                'prof_idx': np.arange(n), 'descr': descr, 'lat': lat, 'lon': lon,
                'YYYYMMDD': ymd, 'year': (ymd // 10000).astype(np.int16),
                'month': ((ymd // 100) % 100).astype(np.int8),
                'bin_a': (np.nan_to_num(bin_a, nan=0).astype(np.int16) if bin_a is not None
                          else np.zeros(n, np.int16))}
        used, fills, outl, cost_out = {}, {}, [], {}
        for V in VARS:
            x, c, w = _read(v, f'prof_{V}'), _read(v, f'prof_{V}clim'), _read(v, f'prof_{V}weight')
            sf = _read(v, f'prof_{V}spatial_scaling_factor')
            if x is None or c is None or w is None:
                continue
            if sf is None:
                sf = np.ones(n)
            if gamma is None:
                gamma = np.ones(n)
            fx, fc, fw = _is_fill_value(x), _is_fill_value(c), _is_fill_weight(w)
            fills[V] = {'value_fill': int(fx.sum()), 'clim_fill': int((fc & ~fx).sum()),
                        'weight_fill': int((fw & ~fx & ~fc).sum()),
                        'weight_ge_9000_kept': int(((np.abs(w) >= 9000) & ~fw).sum())}
            ok = ~fx & ~fc & ~fw & (w > 0)
            cost = np.where(ok, (x - c) ** 2 * w, np.nan)
            sf_ok = np.where(np.isfinite(sf) & (sf > 0), sf, 1.0)
            g_ok = np.where(np.isfinite(gamma) & (gamma > 0), gamma, 1.0)
            divisor = {'final': np.ones(n), 'nosf': sf_ok, 'true': sf_ok * g_ok}
            costs = {k: cost / divisor[k][:, None] for k in KINDS}
            cost_out[V] = cost.astype(np.float32)
            cost_out[f'sf{V}'] = sf.astype(np.float32)
            cost_out['gamma'] = gamma.astype(np.float32)

            with np.errstate(invalid='ignore', all='ignore'):
                nused = ok.sum(axis=1)
                prof[f'{V}_n'] = nused.astype(np.int16)
                mean = np.where(nused > 0, np.nansum(cost, axis=1) / np.maximum(nused, 1), np.nan)
                mx = np.where(nused > 0, np.nanmax(np.where(ok, cost, -np.inf), axis=1), np.nan)
                for k in KINDS:
                    prof[f'{V}_mean_{k}'] = mean / divisor[k]
                    prof[f'{V}_max_{k}'] = mx / divisor[k]
                kmax = np.argmax(np.where(ok, cost, -np.inf), axis=1)
                prof[f'{V}_depth_of_max'] = np.where(nused > 0, depth[kmax], np.nan)
            prof[f'sf{V}'] = sf
            prof['gamma'] = gamma

            used[V] = {k: costs[k][ok] for k in KINDS}
            # outlier candidates: largest true costs in this file
            flat = np.where(ok, costs['true'], -np.inf).ravel()
            k = min(N_OUTLIERS, int(ok.sum()))
            if k:
                top = np.argpartition(flat, -k)[-k:]
                for t in top:
                    i, j = divmod(int(t), cost.shape[1])
                    outl.append({'var': V, **{f'cost_{k}': float(costs[k][i, j]) for k in KINDS},
                                 'source': source, 'subtype': subtype, 'file': path.name, 'prof_idx': i,
                                 'descr': descr[i], 'lat': lat[i], 'lon': lon[i], 'YYYYMMDD': int(ymd[i]),
                                 'depth': float(depth[j]), 'value': float(x[i, j]), 'clim': float(c[i, j]),
                                 'diff': float(x[i, j] - c[i, j]), 'weight_final': float(w[i, j]),
                                 'sf': float(sf_ok[i]), 'gamma': float(g_ok[i]),
                                 'weight_true': float(w[i, j] / (sf_ok[i] * g_ok[i]))})

    _write_cost_nc(out_dir / 'costs' / source / f'{path.stem}__cost.nc', path, source, subtype,
                   n, depth, descr, lat, lon, ymd, hms, cost_out)
    return pd.DataFrame(prof), used, fills, outl


def _write_cost_nc(path, src_path, source, subtype, n, depth, descr, lat, lon, ymd, hms, cost_out):
    path.parent.mkdir(parents=True, exist_ok=True)
    with netCDF4.Dataset(path, 'w') as nc:
        nc.setncatts({'description': 'cost = (prof_V - prof_Vclim)^2 * prof_Vweight per level, recomputed '
                                     'from final chain output. NaN = excluded (fill in V, Vclim or Vweight, '
                                     'or Vweight == 0). Weights are final (include the area '
                                     'factor gamma and the step 11 factor sfV). Other versions: '
                                     'cost_V / sfV = gamma*(dV/sigma)^2 (step 11 removed); '
                                     'cost_V / (gamma * sfV) = (dV/sigma)^2 (both removed).',
                      'source_file': str(src_path), 'source': source, 'subtype': subtype,
                      'created': datetime.datetime.now().isoformat(timespec='seconds')})
        nc.createDimension('iPROF', n)
        nc.createDimension('iDEPTH', len(depth))
        def var(name, data, dims, dtype, **attrs):
            """Create variable ``name`` in the output file with ``data`` and attributes (compressed)."""
            x = nc.createVariable(name, dtype, dims, zlib=True, complevel=4)
            x[:] = data
            x.setncatts(attrs)
        var('prof_depth', depth, ('iDEPTH',), 'f4', units='m')
        var('prof_lat', lat, ('iPROF',), 'f8')
        var('prof_lon', lon, ('iPROF',), 'f8')
        var('prof_YYYYMMDD', ymd, ('iPROF',), 'i4')
        var('prof_HHMMSS', hms, ('iPROF',), 'i4')
        d = nc.createVariable('prof_descr', str, ('iPROF',))
        d[:] = np.asarray(descr, dtype=object)
        for V in VARS:
            if V in cost_out:
                var(f'cost_{V}', cost_out[V], ('iPROF', 'iDEPTH'), 'f4',
                    long_name=f'({V} - {V}clim)^2 * {V}weight (final weights)')
                var(f'sf{V}', cost_out[f'sf{V}'], ('iPROF',), 'f4',
                    long_name=f'step 11 spatial scaling factor for {V}')
        if 'gamma' in cost_out:
            var('gamma', cost_out['gamma'], ('iPROF',), 'f4',
                long_name='grid-cell area factor prof_area_gamma (applied to T and S weights)')


def _stats_row(arrs):
    """{kind: array} -> {'final_n': n, 'final_min': ..., 'nosf_min': ..., ...}"""
    row = {}
    for k in KINDS:
        for stat, x in _summ(arrs[k]).items():
            if stat == 'n':
                if k == 'final':
                    row['final_n'] = x
            else:
                row[f'{k}_{stat}'] = x
    return row


def _process_group(args):
    """One (source, year) group: all its files."""
    key, files, run_dir, out_dir = args
    profs, outl, per_file = [], [], []
    used = defaultdict(lambda: {k: [] for k in KINDS})
    fills = defaultdict(lambda: defaultdict(int))
    for f in files:
        try:
            df, u, fl, ol = _process_file(f, run_dir, out_dir)
        except Exception as e:
            return key, None, None, None, None, None, f'{f.name}: {type(e).__name__}: {e}'
        profs.append(df)
        for V, arrs in u.items():
            for k in KINDS:
                used[V][k].append(arrs[k])
            per_file.append({'file': f.name, 'source': key[0], 'subtype': df['subtype'].iloc[0] if len(df) else '',
                             'year': key[1], 'var': V, **_stats_row(arrs)})
        for V, d in fl.items():
            for k, x in d.items():
                fills[V][k] += x
        outl += ol
    stats, sample, exact = [], {}, {}
    rng = np.random.default_rng(0)
    for V, lists in used.items():
        arrs = {k: (np.concatenate(lists[k]) if lists[k] else np.array([])) for k in KINDS}
        stats.append({'source': key[0], 'year': key[1], 'var': V, **_stats_row(arrs)})
        # thinned sample (same points for every version) for all-years medians / plots
        n = arrs['final'].size
        idx = np.arange(n) if n <= 50_000 else rng.choice(n, 50_000, replace=False)
        sample[V] = {k: arrs[k][idx] for k in KINDS}
        exact[V] = {k: (arrs[k].size, arrs[k].sum(), arrs[k].min() if n else np.nan,
                        arrs[k].max() if n else np.nan) for k in KINDS}
    outl = sorted(outl, key=lambda r: -r['cost_true'])[:N_OUTLIERS]
    return key, pd.concat(profs, ignore_index=True), stats, per_file, \
        {V: dict(d) for V, d in fills.items()}, (outl, sample, exact), None


def build(run_dir=None, out_dir=None, n_workers=None, sources=None):
    """Compute the per-point costs of a run and store them for the report.

    For every yearly file and variable (T, S), three cost versions per point: final (as the model
    sees it), without the step 11 factor, and the "true" cost (ΔV/σ)², with fills excluded. Writes to
    ``out_dir``: ``profiles.parquet`` (per-profile summary), ``stats_source_year.csv``,
    ``stats_file.csv``, ``stats_source.csv``, ``fills.csv``, ``outliers.csv``, ``meta.json`` and
    per-file NetCDF cost files under ``costs/<source>/``.

    Parameters
    ----------
    run_dir : str or Path, optional
        The run folder (default: the newest finished step 11 run).
    out_dir : str or Path, optional
        Output folder (default: ``<diagnostics.cost_report_dir>/cost_<run stamp>__<now>``).
    n_workers : int, optional
        Parallel workers (default: min(4, CPUs - 1)).
    sources : list of str, optional
        Only these sources (for quick tests).

    Returns
    -------
    Path
        The output folder.
    """
    run_dir = common.yearly_dir(run_dir) if run_dir else common.latest_processed_run()
    stamp = datetime.datetime.now().strftime('%Y%m%d_%H%M%S')
    out_dir = Path(out_dir) if out_dir else REPORT_PARENT / f'cost_{common.run_timestamp(run_dir)}__{stamp}'
    out_dir.mkdir(parents=True, exist_ok=True)
    n_workers = n_workers or min(4, max(1, (os.cpu_count() or 2) - 1))
    print(f'Run: {run_dir}\nOut: {out_dir}')

    groups = defaultdict(list)
    for f in sorted(run_dir.rglob('*.nc')):
        src, _ = common.classify(f, run_dir)
        if sources and src not in sources:
            continue
        groups[(src, common.file_year(f))].append(f)
    print(f'{sum(map(len, groups.values()))} files in {len(groups)} (source, year) groups, {n_workers} workers')

    profs, stats, per_file, outl = [], [], [], []
    fills = defaultdict(lambda: defaultdict(lambda: defaultdict(int)))
    samples = defaultdict(lambda: defaultdict(lambda: {k: [] for k in KINDS}))
    exact = defaultdict(lambda: defaultdict(list))
    errors = []
    t0 = time.time()
    # largest groups first so the pool stays busy
    order = sorted(groups.items(), key=lambda kv: -sum(f.stat().st_size for f in kv[1]))
    with concurrent.futures.ProcessPoolExecutor(max_workers=n_workers) as ex:
        futs = [ex.submit(_process_group, (k, fs, run_dir, out_dir)) for k, fs in order]
        for i, fut in enumerate(concurrent.futures.as_completed(futs), 1):
            key, df, st, pf, fl, extra, err = fut.result()
            if err:
                errors.append(err); print(f'  [error] {err}')
                continue
            profs.append(df); stats += st; per_file += pf
            ol, sample, ex_ = extra
            outl += ol
            for V, d in fl.items():
                for k, x in d.items():
                    fills[key[0]][V][k] += x
            for V, d in sample.items():
                for k in KINDS:
                    samples[key[0]][V][k].append(d[k])
            for V, e in ex_.items():
                exact[key[0]][V].append(e)
            print(f'  [{i}/{len(futs)}] {key[0]} {key[1]}  ({time.time() - t0:.0f}s)', flush=True)

    pd.concat(profs, ignore_index=True).to_parquet(out_dir / 'profiles.parquet', index=False)
    pd.DataFrame(stats).sort_values(['var', 'source', 'year']).to_csv(out_dir / 'stats_source_year.csv', index=False)
    pd.DataFrame(per_file).sort_values(['var', 'source', 'year', 'file']).to_csv(out_dir / 'stats_file.csv', index=False)

    # all-years per source: exact n/min/max/mean; median from the per-group samples
    rows = []
    for src in sorted(exact):
        for V in VARS:
            if V not in exact[src]:
                continue
            groups_e = exact[src][V]
            n = sum(e['final'][0] for e in groups_e)
            row = {'source': src, 'var': V, 'final_n': int(n)}
            if n:                          # e.g. XBT S: variable present but never used
                for k in KINDS:
                    e = np.array([g[k] for g in groups_e], dtype=float)
                    smp = np.concatenate(samples[src][V][k])
                    row.update({f'{k}_min': np.nanmin(e[:, 2]), f'{k}_max': np.nanmax(e[:, 3]),
                                f'{k}_mean': e[:, 1].sum() / n,
                                f'{k}_median~': float(np.median(smp)) if smp.size else np.nan})
            rows.append(row)
    pd.DataFrame(rows).to_csv(out_dir / 'stats_source.csv', index=False)
    np.savez_compressed(out_dir / 'samples.npz', **{f'{s_}__{V}__{k}': np.concatenate(samples[s_][V][k])
                                                   for s_ in samples for V in samples[s_] for k in KINDS})

    fr = [{'source': s, 'var': V, **d} for s in sorted(fills) for V, d in fills[s].items()]
    pd.DataFrame(fr).to_csv(out_dir / 'fills.csv', index=False)
    od = pd.DataFrame(outl)
    if len(od):
        od = od.sort_values('cost_true', ascending=False).groupby('var').head(N_OUTLIERS)
    od.to_csv(out_dir / 'outliers.csv', index=False)
    meta = {'run': str(run_dir), 'created': stamp, 'errors': errors,
            'n_files': sum(map(len, groups.values())), 'sources_filter': sources}
    (out_dir / 'meta.json').write_text(json.dumps(meta, indent=2))
    print(f'Done in {time.time() - t0:.0f}s -> {out_dir}')
    return out_dir


if __name__ == '__main__':
    ap = argparse.ArgumentParser(description='Recompute T/S cost vs climatology from final chain files.')
    ap.add_argument('--processed', default=None, help='Run dir (default: newest finished step-11 run)')
    ap.add_argument('--out', default=None)
    ap.add_argument('-n', '--n_workers', type=int, default=None)
    ap.add_argument('--sources', nargs='*', default=None, help='Only these sources (quick tests)')
    a = ap.parse_args()
    build(a.processed, a.out, a.n_workers, a.sources)
