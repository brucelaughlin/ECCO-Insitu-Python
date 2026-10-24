#!/usr/bin/env python
"""
Rewrite MITgcm pkg/profiles input netCDF files with compact variable types
and netCDF-4 (HDF5) zlib compression, with sanity checks before and after.

Steps for each input file
-------------------------
  1. Scan every variable (counts, fill values, min/max/mean, integer-ness).
  2. Print the type plan: input type -> output type for every variable.
  3. Pre-conversion sanity check (table of ERRORs and WARNings).
     Any ERROR aborts the conversion of that file; nothing is written.
  4. Convert (written to <output>.partial, renamed when verified).
  5. Post-conversion check: every value is read back and compared with the
     source, statistics of all fields are tabulated, and the legality
     checks are re-run on the output.

Type rules
----------
  prof_YYYYMMDD, prof_HHMMSS, prof_bin_id_a  -> int32   (read by MITgcm)
  prof_point, prof_basin, prof_bin_id_*,
  prof_*flag, prof_*_code                    -> int32
  prof_interp_*                              -> unchanged (kept float64)
  prof_date                                  -> float64, recomputed as the
      MATLAB datenum of prof_YYYYMMDD + prof_HHMMSS (MITprof files have been
      found with the month number stored in its minutes)
  any other float64 variable                 -> float32
      unless float32 would be unsafe for it: values beyond float32 range,
      or integer-valued data larger than 2**24 (then kept float64)
  integer / char variables                   -> unchanged
      (int64 is narrowed to int32 when the values fit)

Sanity checks (ERROR aborts the conversion, WARN is reported only)
------------------------------------------------------------------
  ERROR  int32 variables and prof_interp_i/j: NaN/Inf, fractional values,
         or values outside int32 range
  ERROR  prof_YYYYMMDD, prof_HHMMSS, prof_lon, prof_lat: NaN/Inf or fill
         values (every profile needs a legal date, time and position)
  ERROR  prof_YYYYMMDD not a valid calendar date or outside
         19000101-21001231; prof_HHMMSS not a valid 24-hour clock time
         (000000-235959, minutes and seconds 00-59)
  ERROR  prof_lat outside [-90, 90]; prof_lon outside [-180, 360];
         NaN/Inf in lon, lat, depth, weights or prof_interp_* variables
  ERROR  negative depth; negative weight (other than the fill value)
  ERROR  weight > 0 where the matching data value is NaN or fill
  ERROR  prof_T / prof_S without a matching weight variable
  ERROR  prof_interp_weights outside [0, 1] or not summing to 1 +- 1e-4
         (the test MITgcm applies); negative prof_interp_i/j
  ERROR  (output only) prof_date differs from prof_YYYYMMDD/prof_HHMMSS
         by more than 0.5 s; on input this is a WARN, since it is recomputed
  WARN   depth not increasing; weights equal to the fill value; weighted
         prof_T outside [-3, 40] or prof_S outside [0, 45]; NaN/Inf in
         other variables; tiny values that underflow to 0 in float32

pkg/profiles reads everything with NF_GET_VARA_DOUBLE, so the types are
converted back to double by the netCDF library at read time. No MITgcm code
change is needed, but the netCDF library that MITgcm links against must
support netCDF-4/HDF5 (check with `nc-config --has-nc4`).

Variables dimensioned by iPROF are chunked as (chunk_prof, <all other dims>)
to match pkg/profiles, which reads blocks of up to 1000 whole profiles
(profiles_readvector.F).

Requires numpy and netCDF4 (conda install -c conda-forge netcdf4).

Examples
--------
  python profiles_compact.py argo_2010.nc                 # -> argo_2010_compact.nc
  python profiles_compact.py argo_2010.nc --check-only    # plan + checks only
  python profiles_compact.py *.nc --suffix _z --complevel 4
  python profiles_compact.py in.nc -o out.nc --chunk-prof 500
"""

import argparse
import datetime
import glob
import os
import re
import sys

import numpy as np
import netCDF4

# integers that MITgcm reads and truncates with INT()
INT_VARS = ('prof_YYYYMMDD', 'prof_HHMMSS', 'prof_bin_id_a')
# other integer quantities written by the MITprof processing chain
INT_PATTERN = re.compile(r'^prof_(point|basin|bin_id_\w+|\w*flag|\w+_code)$')
# integer-valued (MITgcm applies INT()) but kept float64 like all prof_interp_*
INTERP_INT_VARS = ('prof_interp_i', 'prof_interp_j')
KEEP_PREFIXES = ('prof_interp',)
# MATLAB datenum (fractional days, ~7e5; float32 would give ~1.5 hour
# resolution). Not read by MITgcm. Recomputed from prof_YYYYMMDD and
# prof_HHMMSS because MITprof files have been found where its minutes hold
# the month number (datestr 'mm' vs 'MM' mix-up).
DATE_VAR = 'prof_date'
DATENUM_1970 = 719529          # MATLAB datenum(1970, 1, 1)
DATE_TOL_SEC = 0.5             # prof_date vs YYYYMMDD/HHMMSS tolerance
WEIGHT_PATTERN = re.compile(r'^prof_\w+weight$')
# legal range of prof_YYYYMMDD (ERROR outside)
DATE_MIN, DATE_MAX = 19000101, 21001231
# plausible ranges for weighted observations (WARN only)
PHYS_RANGE = {'prof_T': (-3.0, 40.0), 'prof_S': (0.0, 45.0)}
PROF_DIM = 'iPROF'

INT_FILL_DEFAULT = -9999
F32_MAX = float(np.finfo(np.float32).max)
F32_EXACT_INT = 2**24
I32_MIN, I32_MAX = -2**31, 2**31 - 1

# attributes whose type must follow the variable's (new) type
TYPED_ATTRS = ('missing_value', 'valid_min', 'valid_max', 'valid_range',
               'actual_range')
CLASSIC_DTYPES = {np.dtype(t) for t in
                  ('i1', 'i2', 'i4', 'f4', 'f8', 'S1')}

# number of elements read/written per block when streaming a variable
BLOCK_ELEMS = 50_000_000


# ----------------------------------------------------------------------------
# helpers
# ----------------------------------------------------------------------------

def is_int_var(name):
    return name in INT_VARS or bool(INT_PATTERN.match(name))


def open_nc(path, mode='r', **kw):
    ds = netCDF4.Dataset(path, mode, **kw)
    ds.set_auto_maskandscale(False)
    ds.set_auto_chartostring(False)
    return ds


def row_blocks(var, chunk_prof):
    """Yield index keys covering `var` in blocks along its first dim."""
    if var.ndim == 0:
        yield Ellipsis
        return
    n0 = var.shape[0]
    row = int(np.prod(var.shape[1:], dtype=np.int64)) or 1
    step = max(chunk_prof, (BLOCK_ELEMS // row) // chunk_prof * chunk_prof)
    for s in range(0, n0, step):
        yield slice(s, min(s + step, n0))


def fill_values(var):
    """Fill / missing values declared on a variable."""
    vals = []
    for att in ('_FillValue', 'missing_value'):
        if att in var.ncattrs():
            vals.extend(np.atleast_1d(var.getncattr(att)).tolist())
    return vals


def valid_mask(a, fills):
    ok = np.isfinite(a)
    for f in fills:
        ok &= (a != f)
    return ok


def datenum(yyyymmdd, hhmmss):
    """MATLAB datenum of legal prof_YYYYMMDD / prof_HHMMSS values."""
    t = np.asarray(yyyymmdd, dtype=np.float64).astype(np.int64).ravel()
    s = np.asarray(hhmmss, dtype=np.float64).astype(np.int64).ravel()
    y, m, d = t // 10000, (t // 100) % 100, t % 100
    day = (((y - 1970).astype('M8[Y]').astype('M8[M]') + (m - 1))
           .astype('M8[D]') + (d - 1))
    secs = s // 10000 * 3600 + (s // 100) % 100 * 60 + s % 100
    return day.astype(np.int64) + DATENUM_1970 + secs / 86400.0


def fits_int32(x):
    x = float(x)
    return np.isfinite(x) and x == np.trunc(x) and I32_MIN <= x <= I32_MAX


def fmt_num(x, digits=6):
    if x is None:
        return '-'
    if float(x) == np.trunc(x) and abs(x) < 1e15:
        return str(int(x))
    return f'{x:.{digits}g}'


def examples(vals, n=3):
    return ', '.join(fmt_num(v, 12) for v in list(vals)[:n])


def print_table(title, headers, rows, right=()):
    print(f'\n  {title}')
    if not rows:
        print('    (none)')
        return
    w = [max(len(h), *(len(r[i]) for r in rows)) for i, h in enumerate(headers)]

    def line(cells):
        return '    ' + '  '.join(c.rjust(w[i]) if i in right else c.ljust(w[i])
                                  for i, c in enumerate(cells)).rstrip()
    print(line(headers))
    print('    ' + '  '.join('-' * x for x in w))
    for r in rows:
        print(line(r))


# ----------------------------------------------------------------------------
# 1. statistics
# ----------------------------------------------------------------------------

def scan_variable(var, chunk_prof):
    """One pass over a variable. 'valid' = finite and not a fill value."""
    st = dict(numeric=var.dtype.kind in 'fiu', n=int(var.size),
              n_nonfinite=0, n_fill=0, n_valid=0, vmin=None, vmax=None,
              mean=None, n_nonint=0, nonint_ex=[], n_out_i32=0,
              maxabs_all=0.0, n_underflow=0)
    if not st['numeric']:
        return st
    fills = fill_values(var)
    vsum = 0.0
    for key in row_blocks(var, chunk_prof):
        a = np.asarray(var[key], dtype=np.float64).ravel()
        fin = np.isfinite(a)
        st['n_nonfinite'] += int(a.size - np.count_nonzero(fin))
        a = a[fin]
        if not a.size:
            continue
        st['maxabs_all'] = max(st['maxabs_all'], float(np.abs(a).max()))
        v = a[valid_mask(a, fills)]
        st['n_fill'] += int(a.size - v.size)
        if not v.size:
            continue
        st['n_valid'] += int(v.size)
        lo, hi = float(v.min()), float(v.max())
        st['vmin'] = lo if st['vmin'] is None else min(st['vmin'], lo)
        st['vmax'] = hi if st['vmax'] is None else max(st['vmax'], hi)
        vsum += float(v.sum())
        frac = v != np.trunc(v)
        if frac.any():
            st['n_nonint'] += int(np.count_nonzero(frac))
            st['nonint_ex'] = (st['nonint_ex'] + v[frac][:3].tolist())[:3]
        st['n_out_i32'] += int(np.count_nonzero((v < I32_MIN) | (v > I32_MAX)))
        with np.errstate(under='ignore', over='ignore'):
            st['n_underflow'] += int(np.count_nonzero(
                (v != 0) & (v.astype(np.float32) == 0)))
    if st['n_valid']:
        st['mean'] = vsum / st['n_valid']
    return st


# ----------------------------------------------------------------------------
# 2. type plan
# ----------------------------------------------------------------------------

def plan_variable(name, var, st):
    """Output dtype / fill / conversion mode for one variable."""
    src = var.dtype
    spec = dict(dtype=src, mode='copy', fill=None, remap=[],
                int_fill=None, rule='unchanged', notes=[])
    has_fill = '_FillValue' in var.ncattrs()
    if has_fill:
        spec['fill'] = var.getncattr('_FillValue')
    if 'scale_factor' in var.ncattrs() or 'add_offset' in var.ncattrs():
        spec['notes'].append('has scale_factor/add_offset, which pkg/profiles '
                             'ignores; values copied raw')

    if not st['numeric']:
        spec['rule'] = 'char/string, unchanged'
        return spec
    if name.startswith(KEEP_PREFIXES):
        spec['rule'] = 'kept (prof_interp_*)'
        return spec
    if name == DATE_VAR:
        if var.dimensions != (PROF_DIM,):
            spec['rule'] = f'kept (unexpected dims {var.dimensions})'
            return spec
        spec.update(dtype=np.dtype('f8'), mode='replace',
                    rule='recomputed from prof_YYYYMMDD + prof_HHMMSS')
        if has_fill:
            spec['fill'] = np.float64(spec['fill'])
        return spec

    if is_int_var(name):
        spec['rule'] = ('integer, read by MITgcm' if name in INT_VARS
                        else 'integer quantity')
        if src.kind in 'iu' and src.itemsize <= 4:
            return spec
        spec.update(dtype=np.dtype('i4'), mode='int')
        # fill values that cannot be stored as int32 are remapped
        spec['remap'] = [f for f in fill_values(var) if not fits_int32(f)]
        spec['int_fill'] = INT_FILL_DEFAULT
        if has_fill:
            f = np.asarray(spec['fill']).ravel()[0]
            spec['fill'] = np.int32(f if fits_int32(f) else INT_FILL_DEFAULT)
        return spec

    if src.kind in 'iu':
        if src.itemsize == 8 and st['n_out_i32'] == 0 and all(
                fits_int32(f) for f in fill_values(var)):
            spec.update(dtype=np.dtype('i4'), mode='cast', rule='int64 narrowed')
            if has_fill:
                spec['fill'] = np.int32(spec['fill'])
        return spec

    if src == np.float32:
        spec['rule'] = 'already float32'
        return spec

    if st['maxabs_all'] > F32_MAX:
        spec['rule'] = 'kept: exceeds float32 range'
        return spec
    maxabs_valid = max(abs(st['vmin'] or 0.0), abs(st['vmax'] or 0.0))
    if st['n_nonint'] == 0 and maxabs_valid > F32_EXACT_INT:
        spec['rule'] = 'kept: integers > 2**24'
        return spec
    spec.update(dtype=np.dtype('f4'), mode='cast', rule='float64 -> float32')
    if has_fill:
        spec['fill'] = np.float32(spec['fill'])
    return spec


def convert_block(a, spec, key):
    if spec['mode'] == 'int':
        a = np.asarray(a, dtype=np.float64)
        for f in spec['remap']:
            a = np.where((a == f) | (np.isnan(a) & np.isnan(f)),
                         spec['int_fill'], a)
        return a.astype(np.int32)
    if spec['mode'] == 'cast':
        return np.asarray(a).astype(spec['dtype'])
    if spec['mode'] == 'replace':
        return spec['values'][key]
    return a


def convert_attr(val, spec):
    """Convert a typed attribute (missing_value etc.) to the output type."""
    val = np.asarray(val)
    if spec['dtype'].kind == 'i' and val.dtype.kind == 'f':
        val = np.array([v if fits_int32(v) else spec['int_fill']
                        for v in val.ravel()], dtype=np.float64)
    return val.astype(spec['dtype'])


# ----------------------------------------------------------------------------
# 3./5. legality checks (run on the input and again on the output)
# ----------------------------------------------------------------------------

def check_file(ds, stats, chunk_prof, phase='input'):
    """Return a list of (severity, variable, message). phase is 'input' or
    'output'; an inconsistent prof_date is a WARN on input (it will be
    recomputed) and an ERROR on output."""
    issues = []

    def err(v, m):
        issues.append(('ERROR', v, m))

    def warn(v, m):
        issues.append(('WARN', v, m))

    V = ds.variables
    for d in (PROF_DIM, 'iDEPTH'):
        if d not in ds.dimensions:
            err('-', f'missing dimension {d}')
    for n in ('prof_YYYYMMDD', 'prof_HHMMSS', 'prof_lon', 'prof_lat'):
        if n not in V:
            err(n, 'required variable missing')
    depth_name = 'prof_depth' if 'prof_depth' in V else (
        'depth' if 'depth' in V else None)
    if depth_name is None:
        err('prof_depth', 'required variable missing (prof_depth or depth)')

    # per-variable checks from the statistics
    for name, var in V.items():
        st = stats[name]
        if not st['numeric']:
            continue
        must_int = is_int_var(name) or name in INTERP_INT_VARS
        strict = (must_int or name in ('prof_lon', 'prof_lat', depth_name)
                  or WEIGHT_PATTERN.match(name) or name.startswith('prof_interp'))
        if st['n_nonfinite']:
            (err if strict else warn)(name, f"{st['n_nonfinite']} NaN/Inf values")
        if must_int:
            if st['n_nonint']:
                err(name, f"{st['n_nonint']} non-integer values, integers "
                          f"required (e.g. {examples(st['nonint_ex'])})")
            if st['n_out_i32']:
                err(name, f"{st['n_out_i32']} values outside int32 range")
        if name in INTERP_INT_VARS and st['vmin'] is not None and st['vmin'] < 0:
            err(name, f"negative index (min {fmt_num(st['vmin'])})")
        if WEIGHT_PATTERN.match(name):
            if st['vmin'] is not None and st['vmin'] < 0:
                err(name, f"negative weights (min {fmt_num(st['vmin'])})")
            if st['n_fill']:
                warn(name, f"{st['n_fill']} weights equal the fill value "
                           '(MITgcm skips them because they are not > 0)')
        if st['n_underflow'] and not is_int_var(name):
            warn(name, f"{st['n_underflow']} tiny nonzero values underflow "
                       'to 0 in float32')

    # dates, times, positions (1-D, read whole)
    def valid_1d(n):
        a = np.asarray(V[n][:], dtype=np.float64).ravel()
        return a[valid_mask(a, fill_values(V[n]))]

    if 'prof_YYYYMMDD' in V:
        t = np.trunc(valid_1d('prof_YYYYMMDD')).astype(np.int64)
        y, m, d = t // 10000, (t // 100) % 100, t % 100
        mdays = np.array([31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31])
        leap = (y % 4 == 0) & ((y % 100 != 0) | (y % 400 == 0))
        mlen = mdays[np.clip(m, 1, 12) - 1] + (leap & (m == 2))
        bad = (t < 0) | (m < 1) | (m > 12) | (d < 1) | (d > mlen)
        if bad.any():
            err('prof_YYYYMMDD', f'{np.count_nonzero(bad)} invalid dates '
                                 f'(e.g. {examples(t[bad])})')
        odd = ~bad & ((t < DATE_MIN) | (t > DATE_MAX))
        if odd.any():
            err('prof_YYYYMMDD', f'{np.count_nonzero(odd)} dates outside '
                                 f'{DATE_MIN}-{DATE_MAX} '
                                 f'(e.g. {examples(t[odd])})')
    if 'prof_HHMMSS' in V:
        t = np.trunc(valid_1d('prof_HHMMSS')).astype(np.int64)
        bad = ((t < 0) | (t // 10000 > 23) | ((t // 100) % 100 > 59)
               | (t % 100 > 59))
        if bad.any():
            err('prof_HHMMSS', f'{np.count_nonzero(bad)} invalid times '
                               f'(e.g. {examples(t[bad])})')
    for n, lo, hi in (('prof_lat', -90, 90), ('prof_lon', -180, 360)):
        if n in V:
            a = valid_1d(n)
            bad = (a < lo) | (a > hi)
            if bad.any():
                err(n, f'{np.count_nonzero(bad)} values outside [{lo}, {hi}] '
                       f'(e.g. {examples(a[bad])})')
    # every profile needs a real date, time and position
    for n in ('prof_YYYYMMDD', 'prof_HHMMSS', 'prof_lon', 'prof_lat'):
        if n in V and stats[n]['n_fill']:
            a = np.asarray(V[n][:], dtype=np.float64).ravel()
            isf = np.isfinite(a) & ~valid_mask(a, fill_values(V[n]))
            err(n, f"{stats[n]['n_fill']} fill values, a legal value is "
                   f'required (first at iPROF index '
                   f'{int(np.argmax(isf))}, 0-based)')

    # prof_date must agree with prof_YYYYMMDD / prof_HHMMSS
    time_ok = not any(i[0] == 'ERROR' and i[1] in ('prof_YYYYMMDD',
                      'prof_HHMMSS') for i in issues)
    if (DATE_VAR in V and V[DATE_VAR].dimensions == (PROF_DIM,)
            and 'prof_YYYYMMDD' in V and 'prof_HHMMSS' in V and time_ok):
        ref = datenum(V['prof_YYYYMMDD'][:], V['prof_HHMMSS'][:])
        pdv = np.asarray(V[DATE_VAR][:], dtype=np.float64).ravel()
        dsec = np.abs(pdv - ref) * 86400.0
        bad = ~(dsec <= DATE_TOL_SEC)
        if bad.any():
            msg = (f'{np.count_nonzero(bad)} of {bad.size} values differ from '
                   f'prof_YYYYMMDD/prof_HHMMSS by > {DATE_TOL_SEC} s '
                   f'(max {np.nanmax(np.where(bad, dsec, 0)):.0f} s)')
            # signature of the datestr 'mm' (month) vs 'MM' (minutes) mix-up
            fin = bad & np.isfinite(pdv)
            sec = np.round((pdv[fin] - np.floor(pdv[fin])) * 86400).astype(int)
            month = (np.asarray(V['prof_YYYYMMDD'][:], dtype=np.float64)
                     .ravel()[fin].astype(np.int64) // 100) % 100
            if fin.any() and np.all(sec // 60 % 60 == month):
                msg += '; its minutes equal the month number (mm/MM mix-up)'
            if phase == 'input':
                warn(DATE_VAR, msg + '; will be recomputed')
            else:
                err(DATE_VAR, msg)

    if depth_name:
        a = np.asarray(V[depth_name][:], dtype=np.float64).ravel()
        if np.any(a[np.isfinite(a)] < 0):
            err(depth_name, 'negative depths')
        if a.size > 1 and not np.all(np.diff(a) > 0):
            warn(depth_name, 'depths not strictly increasing')

    # observations vs weights
    for n in ('prof_T', 'prof_S'):
        if n in V and n + 'weight' not in V:
            err(n, f'no {n}weight variable')
    for name in V:
        wname = name + 'weight'
        if (wname not in V or not stats[name]['numeric']
                or WEIGHT_PATTERN.match(name)):
            continue
        dv, wv = V[name], V[wname]
        if dv.shape != wv.shape:
            err(name, f'shape {dv.shape} differs from {wname} {wv.shape}')
            continue
        dfill = fill_values(dv)
        n_bad = n_phys = 0
        first, phys_ex = None, []
        lo, hi = PHYS_RANGE.get(name, (-np.inf, np.inf))
        for key in row_blocks(dv, chunk_prof):
            d = np.asarray(dv[key], dtype=np.float64)
            w = np.asarray(wv[key], dtype=np.float64)
            use = w > 0
            ok = valid_mask(d, dfill)
            b = use & ~ok
            if b.any():
                n_bad += int(np.count_nonzero(b))
                if first is None:
                    first = (key.start if isinstance(key, slice) else 0) \
                        + int(np.argwhere(b)[0][0])
            p = use & ok & ((d < lo) | (d > hi))
            if p.any():
                n_phys += int(np.count_nonzero(p))
                phys_ex = (phys_ex + d[p][:3].tolist())[:3]
        if n_bad:
            err(name, f'{n_bad} points with {wname} > 0 but data NaN/fill '
                      f'(first at iPROF index {first}, 0-based)')
        if n_phys:
            warn(name, f'{n_phys} weighted values outside plausible range '
                       f'[{lo:g}, {hi:g}] (e.g. {examples(phys_ex)})')

    # interpolation weights: same tests as profiles_init_fixed.F
    if 'prof_interp_weights' in V:
        st = stats['prof_interp_weights']
        if st['vmin'] is not None and (st['vmin'] < 0 or st['vmax'] > 1):
            err('prof_interp_weights', f"values outside [0, 1] (min "
                f"{fmt_num(st['vmin'])}, max {fmt_num(st['vmax'])})")
        wv = V['prof_interp_weights']
        n_bad, ex = 0, []
        for key in row_blocks(wv, chunk_prof):
            w = np.asarray(wv[key], dtype=np.float64)
            s = w.sum(axis=-1) if w.ndim > 1 else w
            b = ~(np.abs(s - 1.0) <= 1e-4)
            if b.any():
                n_bad += int(np.count_nonzero(b))
                ex = (ex + s[b][:3].tolist())[:3]
        if n_bad:
            err('prof_interp_weights', f'{n_bad} profiles whose weights do '
                f'not sum to 1 +- 1e-4 (e.g. sums {examples(ex)})')
    return issues


def print_issues(title, issues):
    order = {'ERROR': 0, 'WARN': 1}
    rows = [[s, v, m] for s, v, m in sorted(issues, key=lambda i: order[i[0]])]
    print_table(title, ['level', 'variable', 'issue'], rows)
    n_err = sum(1 for i in issues if i[0] == 'ERROR')
    print(f'    -> {n_err} error(s), {len(issues) - n_err} warning(s)')
    return n_err


# ----------------------------------------------------------------------------
# 4./5. conversion and verification
# ----------------------------------------------------------------------------

def write_output(src, dst_path, specs, complevel, chunk_prof, shuffle, fmt,
                 src_path):
    dst = open_nc(dst_path, 'w', format=fmt)
    gattrs = {a: src.getncattr(a) for a in src.ncattrs()}
    stamp = (f'{datetime.date.today().isoformat()}: {os.path.basename(src_path)}'
             f' converted by profiles_compact.py (compact types, zlib '
             f'level {complevel}, shuffle={shuffle}, iPROF chunk {chunk_prof})')
    if any(s['mode'] == 'replace' for s in specs.values()):
        stamp += f'; {DATE_VAR} recomputed from prof_YYYYMMDD and prof_HHMMSS'
    gattrs['history'] = (str(gattrs['history']) + '\n' + stamp
                         if 'history' in gattrs else stamp)
    dst.setncatts(gattrs)

    for dname, dim in src.dimensions.items():
        dst.createDimension(dname, None if dim.isunlimited() else len(dim))

    for name, var in src.variables.items():
        spec = specs[name]
        kw = {}
        if var.ndim > 0 and var.size > 0:
            kw.update(zlib=complevel > 0, complevel=complevel,
                      shuffle=shuffle and complevel > 0)
            if var.dimensions[0] == PROF_DIM:
                kw['chunksizes'] = ((min(chunk_prof, var.shape[0]),)
                                    + tuple(var.shape[1:]))
        if spec['fill'] is not None:
            kw['fill_value'] = spec['fill']
        out = dst.createVariable(name, spec['dtype'], var.dimensions, **kw)
        for att in var.ncattrs():
            if att == '_FillValue':
                continue
            val = var.getncattr(att)
            if att in TYPED_ATTRS and spec['dtype'] != var.dtype:
                val = convert_attr(val, spec)
            out.setncattr(att, val)
        if spec['mode'] == 'replace':
            note = ('recomputed from prof_YYYYMMDD and prof_HHMMSS by '
                    'profiles_compact.py')
            old = var.getncattr('comment') if 'comment' in var.ncattrs() else ''
            out.setncattr('comment', f'{old}; {note}' if old else note)
        for key in row_blocks(var, chunk_prof):
            out[key] = convert_block(var[key], spec, key)
    dst.close()


def verify(src, dst, specs, in_stats, chunk_prof):
    """Compare output with source; return (table rows, output stats, ok)."""
    out_stats, rows, all_ok = {}, [], True
    for name, var in src.variables.items():
        spec, out = specs[name], dst.variables[name]
        st_in = in_stats[name]
        st = out_stats[name] = scan_variable(out, chunk_prof)
        status, dabs, drel = 'OK', 0.0, 0.0
        if out.dtype != spec['dtype'] or out.shape != var.shape:
            status = 'FAIL: dtype/shape'
        else:
            fills = fill_values(var)
            for key in row_blocks(var, chunk_prof):
                a, b = np.asarray(var[key]), np.asarray(out[key])
                expect = convert_block(a, spec, key)
                if not np.array_equal(b, expect,
                                      equal_nan=expect.dtype.kind == 'f'):
                    status = 'FAIL: values differ'
                    break
                if st['numeric']:
                    a64 = a.astype(np.float64)
                    ok = valid_mask(a64, fills)
                    if ok.any():
                        d = np.abs(b.astype(np.float64)[ok] - a64[ok])
                        dabs = max(dabs, float(d.max()))
                        nz = a64[ok] != 0
                        if nz.any():
                            drel = max(drel, float(
                                (d[nz] / np.abs(a64[ok][nz])).max()))
        if status == 'OK' and spec['mode'] == 'replace':
            status = 'OK (recomputed)'
        elif status == 'OK' and st['numeric']:
            for k in ('n_valid', 'n_fill', 'n_nonfinite'):
                if st[k] != st_in[k]:
                    status = f'FAIL: {k} {st_in[k]} -> {st[k]}'
                    break
        all_ok &= status.startswith('OK')
        if st['numeric']:
            rows.append([name, str(out.dtype), str(st['n_valid']),
                         str(st['n_fill']), str(st['n_nonfinite']),
                         fmt_num(st['vmin']), fmt_num(st['vmax']),
                         fmt_num(st['mean']), f'{dabs:.2g}', f'{drel:.2g}',
                         status])
        else:
            rows.append([name, str(out.dtype), str(st['n']), '-', '-', '-',
                         '-', '-', '-', '-', status])
    return rows, out_stats, all_ok


def process_file(src_path, dst_path, args):
    print('=' * 78)
    print(f'{src_path}')
    try:
        src = open_nc(src_path)
    except OSError as e:
        print(f'  ERROR: cannot open: {e}')
        return False
    dims = ', '.join(f'{n}={len(d)}' for n, d in src.dimensions.items())
    print(f'  format {src.file_format}; dimensions: {dims}')

    in_stats = {n: scan_variable(v, args.chunk_prof)
                for n, v in src.variables.items()}
    specs = {n: plan_variable(n, v, in_stats[n])
             for n, v in src.variables.items()}

    rows = []
    for n, v in src.variables.items():
        s = specs[n]
        rule = s['rule'] + ('; ' + '; '.join(s['notes']) if s['notes'] else '')
        rows.append([n, ' x '.join(v.dimensions) or '(scalar)',
                     str(v.dtype), str(s['dtype']), rule])
    print_table('TYPE PLAN', ['variable', 'dims', 'input', 'output', 'rule'],
                rows)

    n_err = print_issues('PRE-CONVERSION SANITY CHECK',
                         check_file(src, in_stats, args.chunk_prof))
    if n_err:
        src.close()
        print(f'\n  ABORTED: {src_path} has {n_err} error(s); not converted')
        return False
    if args.check_only:
        src.close()
        return True

    # dates and times are legal at this point
    spec = specs.get(DATE_VAR)
    if spec is not None and spec['mode'] == 'replace':
        spec['values'] = datenum(src['prof_YYYYMMDD'][:], src['prof_HHMMSS'][:])

    fmt = args.format
    if fmt == 'auto':
        fmt = 'NETCDF4_CLASSIC'
        if any(s['dtype'] not in CLASSIC_DTYPES for s in specs.values()):
            fmt = 'NETCDF4'
    tmp_path = dst_path + '.partial'
    print(f'\n  writing {dst_path} ({fmt}) ...')
    write_output(src, tmp_path, specs, args.complevel, args.chunk_prof,
                 not args.no_shuffle, fmt, src_path)

    dst = open_nc(tmp_path)
    rows, out_stats, ok = verify(src, dst, specs, in_stats, args.chunk_prof)
    print_table('POST-CONVERSION CHECK (output statistics over valid values; '
                'delta = output - input)',
                ['variable', 'type', 'valid', 'fill', 'NaN/Inf', 'min', 'max',
                 'mean', 'max|delta|', 'max rel', 'status'],
                rows, right=(2, 3, 4, 5, 6, 7, 8, 9))
    n_err_out = print_issues('POST-CONVERSION SANITY CHECK (output file)',
                             check_file(dst, out_stats, args.chunk_prof,
                                        phase='output'))
    dst.close()
    src.close()

    if not ok or n_err_out:
        os.remove(tmp_path)
        print(f'\n  FAILED: output did not verify; {tmp_path} removed')
        return False
    os.replace(tmp_path, dst_path)
    s_in, s_out = os.path.getsize(src_path), os.path.getsize(dst_path)
    print(f'\n  OK: {dst_path}\n  size {s_in/2**20:.2f} MiB -> '
          f'{s_out/2**20:.2f} MiB ({100*s_out/s_in:.1f}%)')
    return True


def expand_inputs(args, suffix):
    """Expand wildcards ourselves: PowerShell and cmd.exe pass them through
    unexpanded. Files matched by a wildcard that already carry the output
    suffix (outputs of an earlier run) are skipped."""
    inputs = []
    for a in args:
        if not glob.has_magic(a):
            inputs.append(a)
            continue
        matches = sorted(glob.glob(a))
        if not matches:
            print(f'WARNING: no files match {a}')
        for m in matches:
            if suffix and os.path.splitext(m)[0].endswith(suffix):
                print(f'skipping {m} (name already ends with {suffix})')
            elif m.endswith('.partial'):
                print(f'skipping {m} (incomplete output of an earlier run)')
            else:
                inputs.append(m)
    return inputs


def main():
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('inputs', nargs='+', help='pkg/profiles netCDF file(s)')
    p.add_argument('-o', '--output', help='output file (single input only)')
    p.add_argument('--suffix', default='_compact',
                   help='suffix added to the input name when -o is not given '
                        '(default: %(default)s)')
    p.add_argument('--check-only', action='store_true',
                   help='print the type plan and sanity check; write nothing')
    p.add_argument('--complevel', type=int, default=2, choices=range(0, 10),
                   metavar='0-9', help='zlib level; 0 disables compression '
                   '(default: %(default)s)')
    p.add_argument('--chunk-prof', type=int, default=1000,
                   help='profiles per chunk along iPROF (default: %(default)s)')
    p.add_argument('--no-shuffle', action='store_true',
                   help='disable the HDF5 shuffle filter')
    p.add_argument('--format', default='auto',
                   choices=('auto', 'NETCDF4_CLASSIC', 'NETCDF4'),
                   help='output format; auto picks NETCDF4_CLASSIC unless a '
                        'type needs NETCDF4 (default: %(default)s)')
    p.add_argument('--overwrite', action='store_true',
                   help='overwrite existing output files')
    args = p.parse_args()

    if args.chunk_prof < 1:
        p.error('--chunk-prof must be >= 1')
    inputs = expand_inputs(args.inputs, None if args.output else args.suffix)
    if not inputs:
        p.error('no input files')
    if args.output and len(inputs) > 1:
        p.error('-o/--output can only be used with a single input file')

    failed = []
    for src_path in inputs:
        if args.output:
            dst_path = args.output
        else:
            stem, ext = os.path.splitext(src_path)
            dst_path = stem + args.suffix + (ext or '.nc')
        if not args.check_only:
            if os.path.abspath(dst_path) == os.path.abspath(src_path):
                print(f'ERROR: output would overwrite input {src_path}')
                failed.append(src_path)
                continue
            if os.path.exists(dst_path) and not args.overwrite:
                print(f'ERROR: {dst_path} exists (use --overwrite)')
                failed.append(src_path)
                continue
        if not process_file(src_path, dst_path, args):
            failed.append(src_path)

    if len(inputs) > 1:
        print('=' * 78)
        print(f'{len(inputs) - len(failed)} of {len(inputs)} files OK')
    for f in failed:
        print(f'  FAILED: {f}')
    sys.exit(1 if failed else 0)


if __name__ == '__main__':
    main()
