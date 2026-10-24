#!/usr/bin/env python3
"""
QC checks on NCEI chain output directories.

Four categories of checks:

  1. Internal consistency  — cost outliers, clim coverage, interp displacement,
                             weight-code distribution
  2. Cross-field           — T/S water-mass envelope, T–S weight mismatch,
                             interp position vs declared position
  3. Cross-run regression  — per-source/year key statistics compared between
                             two output directories
  4. External reference    — spot-check prof_Tclim against WOA23 at fixed
                             anchor points

Usage
-----
  # Run all checks on one directory
  ecco-insitu check-qc /path/to/output_dir

  # Also run cross-run regression against a baseline
  ecco-insitu check-qc /path/to/new_dir --baseline /path/to/old_dir

  # Skip the external reference check (no WOA23 file available)
  ecco-insitu check-qc /path/to/output_dir --no-woa

  # Compact: suppress per-file detail, show only summary/warnings
  ecco-insitu check-qc /path/to/output_dir --summary
"""

import argparse
import glob
import os
import sys
import warnings

import numpy as np
import netCDF4

import ecco_insitu.config as config_loader  # noqa: E402

warnings.filterwarnings('ignore', category=RuntimeWarning)

# ---------------------------------------------------------------------------
# constants
# ---------------------------------------------------------------------------

DATENUM_1970 = 719529

# Physical plausibility ranges used for water-mass check (WARN only)
# Derived from WOD global bounds; not meant to be tight
T_RANGE = (-2.5, 40.0)   # degC potential temperature
S_RANGE = (0.0,  42.0)   # PSU

# T–S water-mass envelope: at each depth bin, max salinity for a given T
# (simplified; catches gross inversions and fresh/warm outliers)
# cost thresholds matching the QC step (step06.py since 2026-10-12)
PROFILE_COST_WARN  = 5.0    # mean profile cost — warn above this
PROFILE_COST_ERROR = 16.0   # QC zeroing threshold

# Interp displacement warning threshold (degrees)
INTERP_DISP_WARN = 1.0   # deg; > this suggests grid lookup failure

# Clim spot-check tolerance (degC / PSU)
CLIM_TOL_T = 0.5
CLIM_TOL_S = 0.3

# WOA23 anchor points for external reference check:
# (label, lat, lon, depth_m, month_0indexed, expected_T, expected_S)
WOA23_ANCHORS = [
    ('N.Atl 40N/30W/100m/Aug',  39.5, -30.5, 100, 7,  15.658, None),
    ('Eq.Pac  0N/180W/50m/Jan',  -0.5, 179.5,  50, 0,  27.5,   None),
    ('S.Ocean 55S/0E/200m/Feb', -55.5,   0.5, 200, 1,   2.5,   None),
]

# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------

def open_nc(path):
    ds = netCDF4.Dataset(path, 'r')
    ds.set_auto_maskandscale(False)
    return ds


def valid(a, fill=-9999.0, tol=1e-3):
    """Boolean mask of finite, non-fill values."""
    a = np.asarray(a, dtype=np.float64)
    return np.isfinite(a) & (np.abs(a - fill) > tol)


def source_year(path):
    stem = os.path.basename(path).replace('__ncei_step_10_compact', '').replace('__ncei_step_10', '')
    return os.path.splitext(stem)[0]


def _banner(title, width=72):
    print(f'\n{"="*width}')
    print(f'  {title}')
    print(f'{"="*width}')


def _warn(msg):
    print(f'  WARN  {msg}')


def _error(msg):
    print(f'  ERROR {msg}')


def _ok(msg):
    print(f'  ok    {msg}')


# ---------------------------------------------------------------------------
# 1. Internal consistency checks (per file)
# ---------------------------------------------------------------------------

def check_internal(path, verbose=True):
    """Returns dict of key statistics and issues for one output file."""
    issues = []
    stats  = {}
    ds = open_nc(path)

    nprof = ds.dimensions['iPROF'].size
    stats['nprof'] = nprof

    # --- cost function ---
    # The QC step skips its cost-based zeroing poleward of 60 (dubious_clim_lat_threshold),
    # except for profiles with a high-latitude sigma (prof_Vuncertainty_source
    # 1 = Arctic file, 2 = southern constant), which it does test.  Flag
    # every profile the QC step tested; skip only the ones it exempted.
    HIGH_LAT_THRESH = 60.0
    lat = np.asarray(ds['prof_lat'][:], dtype=np.float64)
    for V in ('T', 'S'):
        if f'prof_{V}cost' not in ds.variables:
            continue
        Vc  = np.asarray(ds[f'prof_{V}cost'][:],   dtype=np.float64)
        Vw  = np.asarray(ds[f'prof_{V}weight'][:], dtype=np.float64)
        use = (Vw > 0) & valid(Vc)
        if not use.any():
            continue
        mean_cost = float(np.nanmean(np.where(use, Vc, np.nan)))
        if Vc.ndim > 1:
            prof_mean_cost = np.where(use.any(axis=1),
                                      np.nanmean(np.where(use, Vc, np.nan), axis=1),
                                      np.nan)
        else:
            prof_mean_cost = np.where(use, Vc, np.nan)
        tested = np.abs(lat) <= HIGH_LAT_THRESH
        src_key = f'prof_{V}uncertainty_source'
        if src_key in ds.variables:
            tested |= np.asarray(ds[src_key][:]) > 0
        n_warn_profs = int(np.nansum((prof_mean_cost > PROFILE_COST_WARN)  & tested))
        n_err_profs  = int(np.nansum((prof_mean_cost > PROFILE_COST_ERROR) & tested))
        stats[f'{V}cost_mean']   = mean_cost
        stats[f'{V}cost_n_warn'] = n_warn_profs
        stats[f'{V}cost_n_err']  = n_err_profs
        if n_err_profs:
            issues.append(f'{n_err_profs}/{nprof} QC-tested profiles with mean {V}cost > {PROFILE_COST_ERROR} (QC threshold missed)')
        elif n_warn_profs:
            issues.append(f'WARN: {n_warn_profs}/{nprof} QC-tested profiles with mean {V}cost > {PROFILE_COST_WARN}')

    # --- climatology coverage ---
    for key in ('prof_Tclim', 'prof_Sclim'):
        if key in ds.variables:
            clim = np.asarray(ds[key][:], dtype=np.float64)
            all_fill = np.all(~valid(clim), axis=1) if clim.ndim > 1 else ~valid(clim)
            n_no_clim = int(all_fill.sum())
            stats[f'{key}_no_clim'] = n_no_clim
            if n_no_clim:
                issues.append(f'WARN: {n_no_clim}/{nprof} profiles with all-fill {key}')

    # --- weight code distribution ---
    for key in ('prof_Tweight_code', 'prof_Sweight_code'):
        if key in ds.variables:
            codes = np.asarray(ds[key][:], dtype=np.int32).ravel()
            unique, counts = np.unique(codes, return_counts=True)
            stats[f'{key}_dist'] = dict(zip(unique.tolist(), counts.tolist()))

    # --- zero-weight profiles ---
    for obs, wt in (('prof_T', 'prof_Tweight'), ('prof_S', 'prof_Sweight')):
        if obs in ds.variables and wt in ds.variables:
            w = np.asarray(ds[wt][:], dtype=np.float64)
            all_zero = np.all(w <= 0, axis=1) if w.ndim > 1 else (w <= 0)
            n_zero = int(all_zero.sum())
            stats[f'{wt}_all_zero'] = n_zero
            if n_zero > nprof * 0.1:
                issues.append(f'WARN: {n_zero}/{nprof} profiles have all-zero {wt}')

    ds.close()
    return stats, issues


# ---------------------------------------------------------------------------
# 2. Cross-field consistency checks (per file)
# ---------------------------------------------------------------------------

def check_cross_field(path, verbose=True):
    issues = []
    ds = open_nc(path)
    nprof = ds.dimensions['iPROF'].size

    # --- T physical range ---
    if 'prof_T' in ds.variables and 'prof_Tweight' in ds.variables:
        T  = np.asarray(ds['prof_T'][:],       dtype=np.float64)
        Tw = np.asarray(ds['prof_Tweight'][:], dtype=np.float64)
        use = (Tw > 0) & valid(T)
        if use.any():
            Tv = T[use]
            n_cold = int((Tv < T_RANGE[0]).sum())
            n_warm = int((Tv > T_RANGE[1]).sum())
            if n_cold:
                issues.append(f'WARN: {n_cold} weighted T values < {T_RANGE[0]}°C '
                               f'(min {Tv.min():.2f})')
            if n_warm:
                issues.append(f'WARN: {n_warm} weighted T values > {T_RANGE[1]}°C '
                               f'(max {Tv.max():.2f})')

    # --- S physical range ---
    if 'prof_S' in ds.variables and 'prof_Sweight' in ds.variables:
        S  = np.asarray(ds['prof_S'][:],       dtype=np.float64)
        Sw = np.asarray(ds['prof_Sweight'][:], dtype=np.float64)
        use = (Sw > 0) & valid(S)
        if use.any():
            Sv = S[use]
            n_fresh = int((Sv < S_RANGE[0]).sum())
            n_salty = int((Sv > S_RANGE[1]).sum())
            if n_fresh:
                issues.append(f'WARN: {n_fresh} weighted S values < {S_RANGE[0]} PSU '
                               f'(min {Sv.min():.3f})')
            if n_salty:
                issues.append(f'WARN: {n_salty} weighted S values > {S_RANGE[1]} PSU '
                               f'(max {Sv.max():.3f})')

    # --- T and S weight consistency ---
    # If both are present, profiles with T weight > 0 but S weight == 0 are expected
    # for some sources; but profiles with S weight > 0 and T weight == 0 are suspicious
    if all(v in ds.variables for v in ('prof_T', 'prof_Tweight', 'prof_S', 'prof_Sweight')):
        Tw = np.asarray(ds['prof_Tweight'][:], dtype=np.float64)
        Sw = np.asarray(ds['prof_Sweight'][:], dtype=np.float64)
        has_T = Tw.max(axis=1) > 0 if Tw.ndim > 1 else Tw > 0
        has_S = Sw.max(axis=1) > 0 if Sw.ndim > 1 else Sw > 0
        s_without_t = int((has_S & ~has_T).sum())
        if s_without_t:
            issues.append(f'WARN: {s_without_t}/{nprof} profiles have S weight but no T weight')

    # --- interp displacement ---
    if all(v in ds.variables for v in ('prof_lon', 'prof_lat',
                                        'prof_interp_lon', 'prof_interp_lat')):
        lon  = np.asarray(ds['prof_lon'][:],        dtype=np.float64)
        lat  = np.asarray(ds['prof_lat'][:],        dtype=np.float64)
        ilon = np.asarray(ds['prof_interp_lon'][:], dtype=np.float64)
        ilat = np.asarray(ds['prof_interp_lat'][:], dtype=np.float64)
        dlon = np.abs(lon - ilon)
        dlon = np.where(dlon > 180, 360 - dlon, dlon)
        dlat = np.abs(lat - ilat)
        n_bad_lon = int((dlon > INTERP_DISP_WARN).sum())
        n_bad_lat = int((dlat > INTERP_DISP_WARN).sum())
        if n_bad_lon:
            issues.append(f'WARN: {n_bad_lon}/{nprof} profiles with interp lon '
                           f'displacement > {INTERP_DISP_WARN}° (max {dlon.max():.2f}°) '
                           f'— possible grid-lookup failure')
        if n_bad_lat:
            issues.append(f'WARN: {n_bad_lat}/{nprof} profiles with interp lat '
                           f'displacement > {INTERP_DISP_WARN}°')

    # --- interp weights sum to 1 ---
    if 'prof_interp_weights' in ds.variables:
        w = np.asarray(ds['prof_interp_weights'][:], dtype=np.float64)
        wsum = w.sum(axis=-1) if w.ndim > 1 else w
        n_bad = int((np.abs(wsum - 1.0) > 1e-4).sum())
        if n_bad:
            issues.append(f'ERROR: {n_bad}/{nprof} profiles with interp weights '
                           f'not summing to 1 ± 1e-4')

    ds.close()
    return issues


# ---------------------------------------------------------------------------
# 3. Cross-run regression (two directories)
# ---------------------------------------------------------------------------

def _file_stats(path):
    """Compact per-file statistics for regression comparison."""
    ds = open_nc(path)
    s = {}
    s['nprof'] = ds.dimensions['iPROF'].size

    for key, wkey in (('prof_T', 'prof_Tweight'), ('prof_S', 'prof_Sweight')):
        if key in ds.variables and wkey in ds.variables:
            v = np.asarray(ds[key][:],  dtype=np.float64)
            w = np.asarray(ds[wkey][:], dtype=np.float64)
            use = (w > 0) & valid(v)
            s[f'{key}_n']    = int(use.sum())
            s[f'{key}_mean'] = float(v[use].mean()) if use.any() else np.nan
        if wkey in ds.variables:
            w = np.asarray(ds[wkey][:], dtype=np.float64)
            all_zero = np.all(w <= 0, axis=1) if w.ndim > 1 else (w <= 0)
            s[f'{wkey}_all_zero'] = int(all_zero.sum())

    for ckey in ('prof_Tclim', 'prof_Sclim'):
        if ckey in ds.variables:
            c = np.asarray(ds[ckey][:], dtype=np.float64)
            vld = valid(c)
            s[f'{ckey}_n']    = int(vld.sum())
            s[f'{ckey}_mean'] = float(c[vld].mean()) if vld.any() else np.nan

    for ckey in ('prof_Tcost', 'prof_Scost'):
        if ckey in ds.variables:
            wkey2 = ckey.replace('cost', 'weight')
            if wkey2 in ds.variables:
                c = np.asarray(ds[ckey][:],  dtype=np.float64)
                w = np.asarray(ds[wkey2][:], dtype=np.float64)
                use = (w > 0) & valid(c)
                s[f'{ckey}_mean'] = float(c[use].mean()) if use.any() else np.nan

    ds.close()
    return s


def check_regression(new_dir, base_dir, summary=False):
    _banner(f'CHECK 3: Cross-run regression\n  new:      {new_dir}\n  baseline: {base_dir}')

    new_files  = {source_year(f): f for f in
                  glob.glob(os.path.join(new_dir,  '**', '*.nc'), recursive=True)}
    base_files = {source_year(f): f for f in
                  glob.glob(os.path.join(base_dir, '**', '*.nc'), recursive=True)}

    common = sorted(set(new_files) & set(base_files))
    only_new  = sorted(set(new_files) - set(base_files))
    only_base = sorted(set(base_files) - set(new_files))

    if only_new:
        print(f'\n  Files only in new ({len(only_new)}):')
        for k in only_new[:10]:
            print(f'    {k}')
        if len(only_new) > 10:
            print(f'    ... and {len(only_new)-10} more')

    if only_base:
        print(f'\n  Files only in baseline ({len(only_base)}):')
        for k in only_base[:10]:
            print(f'    {k}')
        if len(only_base) > 10:
            print(f'    ... and {len(only_base)-10} more')

    # thresholds for flagging regressions
    THRESHOLDS = {
        'nprof':            ('abs', 0),        # must be identical
        'prof_T_n':         ('rel', 0.02),
        'prof_T_mean':      ('abs', 0.1),      # degC
        'prof_S_n':         ('rel', 0.02),
        'prof_S_mean':      ('abs', 0.05),     # PSU
        'prof_Tclim_mean':  ('abs', 0.1),
        'prof_Sclim_mean':  ('abs', 0.05),
        'prof_Tcost_mean':  ('rel', 0.10),
        'prof_Scost_mean':  ('rel', 0.10),
        'prof_Tweight_all_zero': ('abs', 0),
        'prof_Sweight_all_zero': ('abs', 0),
    }

    n_flagged = 0
    for key in common:
        ns = _file_stats(new_files[key])
        bs = _file_stats(base_files[key])
        flags = []
        for stat, (kind, thresh) in THRESHOLDS.items():
            nv = ns.get(stat)
            bv = bs.get(stat)
            if nv is None or bv is None:
                continue
            if np.isnan(nv) or np.isnan(bv):
                continue
            if kind == 'abs':
                diff = abs(nv - bv)
                if diff > thresh:
                    flags.append(f'{stat}: {bv:.4g} → {nv:.4g} (Δ={diff:.4g})')
            else:
                ref = abs(bv) if abs(bv) > 1e-12 else 1e-12
                rdiff = abs(nv - bv) / ref
                if rdiff > thresh:
                    flags.append(f'{stat}: {bv:.4g} → {nv:.4g} ({rdiff*100:.1f}%)')
        if flags:
            n_flagged += 1
            print(f'\n  REGRESSION  {key}')
            for f in flags:
                print(f'    {f}')
        elif not summary:
            print(f'  ok  {key}')

    print(f'\n  {len(common)} files compared, {n_flagged} with regressions, '
          f'{len(only_new)} new-only, {len(only_base)} baseline-only')


# ---------------------------------------------------------------------------
# 4. External reference: spot-check prof_Tclim against WOA23
# ---------------------------------------------------------------------------

def check_woa23_reference(output_dir, woa23_path, summary=False):
    _banner(f'CHECK 4: External reference (WOA23)\n  {woa23_path}')

    woa = open_nc(woa23_path)
    woa_T   = np.asarray(woa['potential_T_monthly'][:], dtype=np.float64)  # (12, depth, lat, lon)
    woa_S   = np.asarray(woa['S_monthly'][:],           dtype=np.float64)
    woa_lat = np.asarray(woa['lat'][:],   dtype=np.float64)
    woa_lon = np.asarray(woa['lon'][:],   dtype=np.float64)
    woa_dep = np.asarray(woa['depth'][:], dtype=np.float64)
    woa.close()

    all_files = sorted(glob.glob(os.path.join(output_dir, '**', '*.nc'), recursive=True))
    if not all_files:
        print('  no output files found')
        return

    n_flagged = 0
    for label, alat, alon, adep, amonth, ref_T, ref_S in WOA23_ANCHORS:
        ilat = int(np.argmin(np.abs(woa_lat - alat)))
        ilon = int(np.argmin(np.abs(woa_lon - alon)))
        idep = int(np.argmin(np.abs(woa_dep - adep)))
        woa_t_val = float(woa_T[amonth, idep, ilat, ilon])
        woa_s_val = float(woa_S[amonth, idep, ilat, ilon])

        # find profiles near this anchor in the output files
        target_month = amonth + 1   # YYYYMMDD months are 1-indexed
        clim_vals = []

        for f in all_files:
            ds = open_nc(f)
            if 'prof_Tclim' not in ds.variables:
                ds.close()
                continue
            lat  = np.asarray(ds['prof_lat'][:],      dtype=np.float64)
            lon  = np.asarray(ds['prof_lon'][:],      dtype=np.float64)
            ymd  = np.asarray(ds['prof_YYYYMMDD'][:], dtype=np.int64)
            clim = np.asarray(ds['prof_Tclim'][:],    dtype=np.float64)
            dep  = np.asarray(ds['prof_depth'][:],    dtype=np.float64) if 'prof_depth' in ds.variables else None
            ds.close()

            month = (ymd % 10000) // 100
            near = ((np.abs(lat - alat) < 2.0) &
                    (np.abs(lon - alon) < 2.0) &
                    (month == target_month))
            if not near.any() or dep is None:
                continue

            # find depth index closest to adep
            idep_obs = int(np.argmin(np.abs(dep - adep)))
            for iprof in np.where(near)[0]:
                v = float(clim[iprof, idep_obs])
                if valid(np.array([v])).any():
                    clim_vals.append(v)
            if len(clim_vals) >= 20:
                break

        if not clim_vals:
            print(f'  {label}: no profiles found near anchor — skipping')
            continue

        clim_mean = float(np.mean(clim_vals))
        diff = abs(clim_mean - woa_t_val)
        status = 'WARN' if diff > CLIM_TOL_T else 'ok  '
        if status == 'WARN':
            n_flagged += 1
        print(f'  {status}  {label}')
        print(f'        WOA23 T={woa_t_val:.3f}°C  '
              f'prof_Tclim mean={clim_mean:.3f}°C  '
              f'|diff|={diff:.3f}°C  '
              f'(n={len(clim_vals)} profiles, tol={CLIM_TOL_T}°C)')

    print(f'\n  {n_flagged}/{len(WOA23_ANCHORS)} anchor points outside tolerance')


# ---------------------------------------------------------------------------
# main driver
# ---------------------------------------------------------------------------

def run_checks_1_2(output_dir, summary=False):
    _banner(f'CHECKS 1+2: Internal consistency & cross-field\n  {output_dir}')

    all_files = sorted(glob.glob(os.path.join(output_dir, '**', '*.nc'), recursive=True))
    if not all_files:
        print('  no output files found')
        return

    total_issues = []
    for path in all_files:
        stats, issues1 = check_internal(path)
        issues2        = check_cross_field(path)
        all_iss = issues1 + issues2
        label = '/'.join(path.split('/')[-2:])
        if all_iss:
            total_issues.extend(all_iss)
            print(f'\n  {label}  ({stats["nprof"]} profiles)')
            for iss in all_iss:
                print(f'    {iss}')
        elif not summary:
            print(f'  ok  {label}  ({stats["nprof"]} profiles)')

    errors = [i for i in total_issues if not i.startswith('WARN')]
    warns  = [i for i in total_issues if i.startswith('WARN')]
    print(f'\n  {len(all_files)} files checked — '
          f'{len(errors)} errors, {len(warns)} warnings')


def cli(argv=None):
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('output_dir', help='NCEI chain output directory to check')
    p.add_argument('--baseline', metavar='DIR',
                   help='prior output directory for cross-run regression (check 3)')
    p.add_argument('--woa23', metavar='FILE',
                   default=None,
                   help='WOA23 full-depth climatology file for check 4 '
                        '(default: paths.climatology_file in the settings)')
    p.add_argument('--no-woa', action='store_true',
                   help='skip check 4 (external WOA23 reference)')
    p.add_argument('--summary', action='store_true',
                   help='suppress per-file OK lines; show only warnings/errors')
    args = p.parse_args(argv)

    if not os.path.isdir(args.output_dir):
        p.error(f'not a directory: {args.output_dir}')
    # a run folder with a yearly/ subfolder: check the yearly files
    for attr in ('output_dir', 'baseline'):
        d = getattr(args, attr, None)
        if d and os.path.isdir(os.path.join(d, 'yearly')):
            setattr(args, attr, os.path.join(d, 'yearly'))

    run_checks_1_2(args.output_dir, summary=args.summary)

    if args.baseline:
        if not os.path.isdir(args.baseline):
            p.error(f'baseline not a directory: {args.baseline}')
        check_regression(args.output_dir, args.baseline, summary=args.summary)

    if args.woa23 is None:
        args.woa23 = config_loader.load()['paths']['climatology_file'] or ''
    if not args.no_woa:
        if not os.path.isfile(args.woa23):
            print(f'\nSkipping check 4: WOA23 file not found: {args.woa23}')
        else:
            check_woa23_reference(args.output_dir, args.woa23, summary=args.summary)

    print()


if __name__ == '__main__':
    sys.exit(cli())
