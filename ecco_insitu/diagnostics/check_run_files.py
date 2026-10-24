#!/usr/bin/env python3
"""
Structural validation of a finished run_ncei.sh output directory.

Complements check_run_qc.py (which checks scientific plausibility) with checks
that the files are complete and internally consistent:

  per file
    - opens, iPROF > 0, every file at the same final step
    - required variables present; prof_interp_* shaped (iPROF, iINTERP=1)
    - weights finite and >= 0; every profile has a nonzero T or S weight
    - no weighted fill: weight > 0 only where the value and its climatology are
      real (not NaN / fill)
    - prof_bin_id_a in 1..10242, prof_bin_id_b in 1..2562 (warning: profiles whose
      nearest LLC90 cell is land get the -9999 fill, i.e. no bin, and are not
      spatially scaled)
    - prof_YYYYMMDD / prof_HHMMSS are valid calendar values and agree with
      prof_date (MATLAB datenum) to within 1 s
    - step 11: prof_Vspatial_scaling_factor present for each V weight, in (0, 1]
  across the run
    - every file on the same depth grid (step 08 trimming off); --allow_mixed_grids
      to accept per-file trimmed grids
  across the run (step 11)
    - for each V, the scaling factors of profiles with V obs sum to exactly 1
      in every (bin, month) — the defining property of 1/n scaling
    - optionally (--input_dir) every input file has an output or a logged
      failure

Usage:
    python diagnostics/check_run_files.py <run_dir> [--input_dir <Interp_Profiles>] [--log <ncei_log>]
Exit code 0 if no errors.
"""

import argparse
import re
import sys
from collections import defaultdict
from pathlib import Path

import netCDF4
import numpy as np

VARS = ('T', 'S')
REQUIRED = ['prof_lon', 'prof_lat', 'prof_date', 'prof_YYYYMMDD', 'prof_HHMMSS', 'prof_depth',
            'prof_bin_id_a', 'prof_bin_id_b']
N_BINS = {'prof_bin_id_a': 10242, 'prof_bin_id_b': 2562}
# Per-profile tile corners stay 1-D by design (tools.MITprof_write_to_nc)
INTERP_1D = {'prof_interp_XC11', 'prof_interp_YC11', 'prof_interp_XCNINJ', 'prof_interp_YCNINJ'}
_DATENUM_1970 = 719529


def _arr(v, i0=None, i1=None):
    x = v[:] if i0 is None else v[i0:i1]
    return np.ma.filled(np.ma.asarray(x, dtype=np.float64), np.nan)


def _datenum(ymd, hms):
    y, m, d = ymd // 10000, (ymd // 100) % 100, ymd % 100
    day = ((((y - 1970).astype('M8[Y]').astype('M8[M]') + (m - 1)).astype('M8[D]') + (d - 1))
           .astype(np.int64))
    return day + _DATENUM_1970 + (hms // 10000 * 3600 + (hms // 100) % 100 * 60 + hms % 100) / 86400


def check_file(path, sums, errors, warnings, grids):
    name = path.name
    err = lambda msg: errors.append(f'{name}: {msg}')
    try:
        nc = netCDF4.Dataset(path)
    except Exception as e:
        err(f'cannot open ({e})')
        return 0
    with nc:
        v = nc.variables
        n = nc.dimensions['iPROF'].size if 'iPROF' in nc.dimensions else 0
        if n == 0:
            err('iPROF == 0')
            return 0
        for r in REQUIRED:
            if r not in v:
                err(f'missing {r}')
        if 'prof_depth' in v:
            grids.setdefault(tuple(np.round(_arr(v['prof_depth']), 3)), []).append(name)
        present = [V for V in VARS if f'prof_{V}weight' in v]
        if not present:
            err('no prof_Tweight / prof_Sweight')
            return n
        for k in v:
            if (k.startswith('prof_interp_') and k not in INTERP_1D
                    and v[k].dimensions != ('iPROF', 'iINTERP')):
                err(f'{k} dims {v[k].dimensions}, expected (iPROF, iINTERP)')
        if 'iINTERP' in nc.dimensions and nc.dimensions['iINTERP'].size != 1:
            err(f"iINTERP = {nc.dimensions['iINTERP'].size}, expected 1")

        any_w = np.zeros(n, bool)
        has = {}
        for V in present:
            w = _arr(v[f'prof_{V}weight'])
            if not np.isfinite(w).all():
                err(f'prof_{V}weight has {int((~np.isfinite(w)).sum())} non-finite values')
            if (w < 0).any():
                err(f'prof_{V}weight has {int((w < 0).sum())} negative values')
            has[V] = (w > 0).any(axis=1) if w.ndim == 2 else (w > 0)
            for what, key in (('value', f'prof_{V}'), ('climatology', f'prof_{V}clim')):
                if key in v:
                    x = _arr(v[key])
                    bad = (w > 0) & (~np.isfinite(x) | (np.abs(x) >= 900))
                    if bad.any():
                        err(f'{int(bad.sum())} points with prof_{V}weight > 0 but fill/NaN {what}')
            any_w |= has[V]
        if (~any_w).any():
            err(f'{int((~any_w).sum())} profiles with no nonzero T or S weight (the zero-weight removal step should drop these)')

        valid_bin = np.ones(n, bool)
        for b, nb in N_BINS.items():
            if b in v:
                ids = _arr(v[b])
                ok = (ids >= 1) & (ids <= nb)
                if b == 'prof_bin_id_a':
                    valid_bin = ok
                bad_val = ~ok & np.isfinite(ids) & (ids != -9999)
                if bad_val.any():
                    err(f'{b}: {int(bad_val.sum())} values outside 1..{nb}')
                elif (~ok).any():
                    warnings.append(f'{name}: {b}: {int((~ok).sum())} profiles with no bin '
                                    f'(nearest LLC90 cell is land); not spatially scaled')

        if all(k in v for k in ('prof_YYYYMMDD', 'prof_HHMMSS', 'prof_date')):
            ymd = _arr(v['prof_YYYYMMDD']).astype(np.int64)
            hms = _arr(v['prof_HHMMSS']).astype(np.int64)
            y, m, d = ymd // 10000, (ymd // 100) % 100, ymd % 100
            hh, mi, ss = hms // 10000, (hms // 100) % 100, hms % 100
            ok_m = (m >= 1) & (m <= 12)
            mdays = np.array([0, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31])
            leap = (y % 4 == 0) & ((y % 100 != 0) | (y % 400 == 0))
            maxd = mdays[np.where(ok_m, m, 1)] + (leap & (m == 2))
            bad = ~(ok_m & (d >= 1) & (d <= maxd) & (hh <= 23) & (mi <= 59) & (ss <= 59) & (hms >= 0))
            if bad.any():
                err(f'{int(bad.sum())} invalid prof_YYYYMMDD/prof_HHMMSS values')
            else:
                diff = np.abs(_datenum(ymd, hms) - _arr(v['prof_date'])) * 86400
                if (diff > 60).any():
                    warnings.append(f'{name}: {int((diff > 60).sum())} profiles where YYYYMMDD/HHMMSS '
                                    f'and prof_date differ by > 60 s (max {diff.max():.0f} s)')

        bins = np.nan_to_num(_arr(v['prof_bin_id_a'])).astype(np.int64) if 'prof_bin_id_a' in v else None
        ym = (_arr(v['prof_YYYYMMDD']).astype(np.int64) // 100) if 'prof_YYYYMMDD' in v else None
        for V in present:
            sname = f'prof_{V}spatial_scaling_factor'
            if sname not in v:
                err(f'missing {sname} (step 11 not applied?)')
                continue
            sf = _arr(v[sname])
            h = has[V] & valid_bin
            if (~np.isfinite(sf[h])).any() or (sf[h] <= 0).any() or (sf[h] > 1 + 1e-12).any():
                err(f'{sname} outside (0, 1] for profiles with {V} obs')
            if bins is not None and ym is not None:
                for key, s in zip(zip(bins[h], ym[h]), sf[h]):
                    sums[V][key] += s
    return n


def main():
    ap = argparse.ArgumentParser(description='Structural validation of an NCEI run directory.')
    ap.add_argument('run_dir')
    ap.add_argument('--input_dir', default=None, help='Check every input file has an output or failure')
    ap.add_argument('--log', default=None, help='run_ncei.sh log, to match missing outputs to failures')
    ap.add_argument('--allow_mixed_grids', action='store_true',
                    help='Do not flag files on different depth grids (step 08 trimming on)')
    a = ap.parse_args()
    run = Path(a.run_dir)
    if (run / 'yearly').is_dir():               # run folder with a yearly/ subfolder
        run = run / 'yearly'
    files = sorted(run.rglob('*.nc'))
    errors, warnings = [], []
    sums = {V: defaultdict(float) for V in VARS}
    grids = {}

    steps = {re.search(r'__ncei_step_(\d+)', f.name).group(1) if '__ncei_step_' in f.name else '?'
             for f in files}
    if len(steps) != 1:
        errors.append(f'files at mixed steps: {sorted(steps)}')
    print(f'{run.name}: {len(files)} files, final step(s) {sorted(steps)}')

    n_prof = 0
    for k, f in enumerate(files, 1):
        n_prof += check_file(f, sums, errors, warnings, grids)
        if k % 50 == 0:
            print(f'  checked {k}/{len(files)}', flush=True)
    print(f'  {n_prof:,} profiles')
    desc = '; '.join(f'{len(fs)} files on {len(g)} levels ({g[0]:g}–{g[-1]:g} m)'
                     for g, fs in sorted(grids.items(), key=lambda kv: -len(kv[1])))
    print(f'  depth grids: {len(grids)} distinct — {desc[:300]}')
    if len(grids) > 1 and not a.allow_mixed_grids:
        errors.append(f'files on {len(grids)} different depth grids (expected one; step 08 trimming on?)')

    for V in VARS:
        if not sums[V]:
            continue
        vals = np.array(list(sums[V].values()))
        bad = np.abs(vals - 1) > 1e-6
        msg = (f'step 11 {V}: {len(vals):,} (bin, month) cells, scaling factors sum to 1 '
               f'in {len(vals) - bad.sum():,}')
        print('  ' + msg)
        if bad.any():
            errors.append(f'{msg}; {int(bad.sum())} cells off (range {vals[bad].min():.4f}–{vals[bad].max():.4f})')

    if a.input_dir:
        inp = Path(a.input_dir)
        made = {re.sub(r'__ncei_step_\d+$', '', f.relative_to(run).with_suffix('').as_posix()) for f in files}
        missing = [f.relative_to(inp).with_suffix('').as_posix() for f in sorted(inp.rglob('*.nc'))
                   if f.relative_to(inp).with_suffix('').as_posix() not in made]
        logged = set()
        if a.log:
            logged = set(re.findall(r'NO OUTPUT FILE WRITTEN for: (\S+)\.nc', Path(a.log).read_text()))
        for m in missing:
            (warnings if Path(m).name in logged else errors).append(
                f'no output for input {m}' + (' (logged chain failure)' if Path(m).name in logged else ''))
        print(f'  inputs: {len(missing)} without output' + (f', {len(logged)} logged failures' if a.log else ''))

    print(f'\nWARNINGS ({len(warnings)}):')
    for w in warnings:
        print('  ' + w)
    print(f'ERRORS ({len(errors)}):')
    for e in errors:
        print('  ' + e)
    print('\nRESULT:', 'PASS' if not errors else 'FAIL')
    sys.exit(1 if errors else 0)


if __name__ == '__main__':
    main()
