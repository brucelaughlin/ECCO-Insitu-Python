"""
Build the per-profile index that every report plot reads.

One row per profile, for two stages:
  raw        — files as delivered by UCSD (Interp_Profiles/)
  processed  — a finished NCEI run (step 11 output)

Columns:
  source, subtype, file      categorical (see common.classify)
  year, month                from prof_YYYYMMDD (file year / month 0 if invalid)
  lat, lon                   float32
  bin_a, bin_b               10242- / 2562-bin ids (0 = unassigned)
  hasT, hasS                 profile carries T / S:
                               raw       — any valid value (non-empty profile,
                                           as on the UCSD page)
                               processed — any nonzero prof_Vweight, i.e. the
                                           profile is actually used for V
                                           (climatology-filled S has zero
                                           weight and is not counted)
  sfT, sfS                   processed only: step 11 spatial scaling factor
                             (NaN if absent; the legacy single factor fills both)

Raw files have no bin ids, so they are binned here with the same
nearest-LLC90-cell lookup as step 02.

Run standalone:
    python -m ecco_insitu.diagnostics.counts_index <report_dir> [--processed RUN] [-n 6]
"""

import argparse
import concurrent.futures
import datetime
import json
import os
import time
from pathlib import Path

import netCDF4
import numpy as np
import pandas as pd

import ecco_insitu.diagnostics.common as common

_CHUNK = 50_000          # profiles read at a time (bounds memory on big files)
_lookup = None           # per-worker BinLookup (raw stage only)


def _worker_init(need_bins):
    global _lookup
    _lookup = common.BinLookup() if need_bins else None


def _read_any(var, i0, i1, test):
    """(i1-i0,) bool: any level of var[i0:i1] passes test."""
    x = var[i0:i1]
    x = np.ma.filled(np.ma.asarray(x, dtype=np.float64), np.nan)
    ok = test(x)
    return ok.any(axis=1) if ok.ndim == 2 else ok


def _valid_T(x):
    return np.isfinite(x) & (x > -99)


def _valid_S(x):
    return np.isfinite(x) & (x > 0)


def _positive(x):
    return np.isfinite(x) & (x > 0)


def _index_file(args):
    path, root, stage = args
    source, subtype = common.classify(path, root)
    try:
        nc = netCDF4.Dataset(path)
    except Exception as e:
        return None, f"{path.name}: cannot open ({e})"
    try:
        n = nc.dimensions['iPROF'].size if 'iPROF' in nc.dimensions else 0
        out = {'source': np.full(n, source, dtype=object),
               'subtype': np.full(n, subtype, dtype=object),
               'file': np.full(n, path.name, dtype=object)}
        if n == 0:
            return pd.DataFrame(out), None

        v = nc.variables
        ymd = np.ma.filled(np.ma.asarray(v['prof_YYYYMMDD'][:], dtype=np.float64), -1).astype(np.int64)
        good = ymd > 0
        out['year']  = np.where(good, ymd // 10000, common.file_year(path)).astype(np.int16)
        out['month'] = np.where(good, (ymd // 100) % 100, 0).astype(np.int8)
        lat = np.ma.filled(np.ma.asarray(v['prof_lat'][:], dtype=np.float64), np.nan)
        lon = np.ma.filled(np.ma.asarray(v['prof_lon'][:], dtype=np.float64), np.nan)
        out['lat'], out['lon'] = lat.astype(np.float32), lon.astype(np.float32)

        if stage == 'raw':
            bins = _lookup.bins(lon, lat)
            for g in common.GRIDS:
                out[f'bin_{g}'] = bins[g].astype(np.int16)
            src = {'T': ('prof_T', _valid_T), 'S': ('prof_S', _valid_S)}
        else:
            for g, (name, _, _) in common.GRIDS.items():
                b = np.ma.filled(np.ma.asarray(v[name][:], dtype=np.float64), 0) if name in v else np.zeros(n)
                out[f'bin_{g}'] = np.nan_to_num(b).astype(np.int16)
            src = {'T': ('prof_Tweight', _positive), 'S': ('prof_Sweight', _positive)}
            legacy = (np.ma.filled(np.ma.asarray(v['prof_spatial_scaling_factor'][:], dtype=np.float64), np.nan)
                      if 'prof_spatial_scaling_factor' in v else None)
            for var in common.VARS:
                name = f'prof_{var}spatial_scaling_factor'
                if name in v:
                    sf = np.ma.filled(np.ma.asarray(v[name][:], dtype=np.float64), np.nan)
                elif legacy is not None:
                    sf = legacy
                else:
                    sf = np.full(n, np.nan)
                out[f'sf{var}'] = sf.astype(np.float32)

        for var, (name, test) in src.items():
            has = np.zeros(n, dtype=bool)
            if name in v:
                for i0 in range(0, n, _CHUNK):
                    i1 = min(n, i0 + _CHUNK)
                    has[i0:i1] = _read_any(v[name], i0, i1, test)
            out[f'has{var}'] = has
        return pd.DataFrame(out), None
    except Exception as e:
        return None, f"{path.name}: {type(e).__name__}: {e}"
    finally:
        nc.close()


def _scan(stage, root, n_workers, sources=None):
    files = sorted(Path(root).rglob('*.nc'))
    if sources:
        files = [f for f in files if f.relative_to(root).parts[0] in sources]
    print(f"\n[{stage}] {len(files)} files under {root}")
    parts, errors = [], []
    t0 = time.time()
    with concurrent.futures.ProcessPoolExecutor(
            max_workers=n_workers, initializer=_worker_init,
            initargs=(stage == 'raw',)) as ex:
        futs = [ex.submit(_index_file, (f, Path(root), stage)) for f in files]
        for k, fut in enumerate(concurrent.futures.as_completed(futs), 1):
            df, err = fut.result()
            if err:
                errors.append(err)
                print(f"  [warn] {err}")
            elif df is not None:
                parts.append(df)
            if k % 25 == 0 or k == len(files):
                print(f"  [{stage}] {k}/{len(files)} files, {time.time() - t0:.0f}s", flush=True)
    df = pd.concat(parts, ignore_index=True) if parts else pd.DataFrame()
    for c in ('source', 'subtype', 'file'):
        if c in df:
            df[c] = df[c].astype('category')
    return df, len(files), errors


def build(report_dir, raw_root=common.RAW_ROOT, processed=None, n_workers=None, sources=None):
    report_dir = Path(report_dir)
    processed = common.yearly_dir(processed) if processed else common.latest_processed_run()
    n_workers = n_workers or min(6, max(1, (os.cpu_count() or 2) - 1))
    raw_p, proc_p, meta_p = common.index_paths(report_dir)
    raw_p.parent.mkdir(parents=True, exist_ok=True)

    meta = {'created': datetime.datetime.now().isoformat(timespec='seconds'),
            'raw_root': str(raw_root), 'processed_run': str(processed),
            'sources_filter': sources}
    for stage, root, path in (('raw', raw_root, raw_p), ('processed', processed, proc_p)):
        df, n_files, errors = _scan(stage, root, n_workers, sources)
        df.to_parquet(path, index=False)
        meta[stage] = {'files': n_files, 'profiles': int(len(df)), 'errors': errors}
        print(f"  [{stage}] {len(df):,} profiles -> {path}")
    meta_p.write_text(json.dumps(meta, indent=2))
    return meta


if __name__ == '__main__':
    ap = argparse.ArgumentParser(description='Build the per-profile index for the in situ report.')
    ap.add_argument('report_dir')
    ap.add_argument('--raw_root', default=str(common.RAW_ROOT))
    ap.add_argument('--processed', default=None,
                    help='Processed run dir (default: newest finished step-11 run)')
    ap.add_argument('-n', '--n_workers', type=int, default=None)
    ap.add_argument('--sources', nargs='*', default=None,
                    help='Only these top-level source dirs (for quick tests)')
    a = ap.parse_args()
    build(a.report_dir, a.raw_root, a.processed, a.n_workers, a.sources)
