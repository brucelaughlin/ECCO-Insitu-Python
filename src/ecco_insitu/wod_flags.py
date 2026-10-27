"""
WOD quality flags for the WOD profiles (``prof_Tflag``, ``prof_Sflag``), applied by step 06.

The input files (``Interp_Profiles``) carry no quality flags, so values that the World Ocean
Database flags as bad would get full weight. This module adds WOD's own per-level flags, as the
original MATLAB chain had them (it read WOD's standard-level exports, whose flags it copied into
``prof_Tflag`` / ``prof_Sflag``; step 06 test 2 then zeroes the weight wherever a flag is > 0).

**Source:** WOD's standard-level files, one per dataset and year, from NCEI::

    https://www.ncei.noaa.gov/data/oceans/woa/WOD/YEARLY/<DATASET>/STD/<DATASET>S<year>.gz

(WOD native format; datasets XBT, CTD, OSD, MRB, GLD, ...; about 2 GB for 1992 on). Their per-level
flags are WOD's standard-deviation checks against its climatology (0 = accepted; 1-9 = rejected, see
WOD's documentation). On the XBT casts of 1992-1998 they agree with the flags in ECCO v4r6's input
files at 99.4% of levels; the rest are WOD's updates since 2018.

**Matching:** by WOD unique cast number, from ``prof_descr`` (``WOD-XBT//<cast>:...``,
``WOD-MRB//<cast>_...``). A profile without a WOD cast number, or not found, gets no flags.

**Depths:** the input files are on their own depth grid, not WOD's standard levels. A depth gets the
flag of the WOD standard level at the same depth, or else the larger of the flags of the two WOD levels
around it (a value interpolated from a flagged level is flagged). Only the per-level flags are used, as
in the original chain; WOD's whole-cast flag is not.

**Depths beyond the cast's WOD levels (DECISION PENDING APPROVAL; setting** ``chain.wod_flags_extend_levels``
**).** Our files can have values deeper than the cast's deepest WOD standard level (or shallower than its
shallowest), where WOD gives no flag. The original chain never had such values (it worked on WOD's
standard levels). Rule: such a value is flagged if the cast's ``N`` end levels nearest to it (the deepest
``N`` for a deeper value, the shallowest ``N`` for a shallower one) are **all** flagged; it then gets the
flag of the nearest level. Otherwise it gets 0.

- ``N = 2`` (default, proposed): a flagged bottom is extended only when it is more than one bad level.
- ``N = 1``: extend whenever the last level is flagged.
- ``N = 0``: never extend (values beyond the cast's WOD levels are unflagged).

Effect on Interp_Profiles (counted 2026-10-27), extra values flagged beyond ``N = 0``:

==========  =================  =================
source      ``N = 1``          ``N = 2``
==========  =================  =================
XBT T       8,452 (0.023%)     7,127 (0.020%)
CTD T / S   8,650 / 12,412     5,628 / 7,499
MRB T / S   6,786 / 5,481      4,201 / 2,231 (S: 0.093%)
GLD T / S   25,877 / 40,116    17,223 / 26,475
==========  =================  =================

For comparison, ``N = 0`` already flags 1.1-3.1% of each source's values. See ``docs/known_issues.md``
(KI-19).

**Setup** (once; see ``paths.wod_std_dir`` and ``chain.apply_wod_flags`` in the settings)::

    ecco-insitu wod-flags --download     # the standard-level files, into <wod_std_dir>/<dataset>/
    ecco-insitu wod-flags --build        # read them once into <wod_std_dir>/index/ (fast lookup)
    ecco-insitu wod-flags                # check

``--build`` reads the files with ``wodpy`` (``pip install wodpy``). A chain run with
``chain.apply_wod_flags: true`` checks the index first and stops with these instructions if it is missing.
"""
import argparse
import datetime
import gzip
import os
import re
import shutil
import sys
import tempfile
import urllib.error
import urllib.request
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np

BASE = 'https://www.ncei.noaa.gov/data/oceans/woa/WOD/YEARLY'
#: WOD datasets searched.
DATASETS = ('XBT', 'CTD', 'OSD', 'MRB', 'GLD')
#: The datasets searched first, by the instrument code in ``prof_descr`` (``WOD-<code>//``); the rest of
#: ``DATASETS`` is searched after them, since a cast can sit in another dataset than its code says (e.g.
#: some 'CTD' casts of 1995 are in WOD's OSD files).
CODE_TO_DATASETS = {'XBT': ('XBT',), 'CTD': ('CTD', 'OSD'), 'OSD': ('OSD', 'CTD'), 'XCTD': ('CTD', 'XBT', 'OSD'),
                    'STD': ('OSD', 'CTD'), 'MRB': ('MRB',), 'GLD': ('GLD',)}
FIRST_YEAR = 1992
_CAST = re.compile(r'^WOD-([A-Z]+)//(\d+)')
#: WOD variable codes: 1 temperature, 2 salinity.
_VAR_CODES = {'prof_T': 1, 'prof_S': 2}


def std_file(wod_std_dir, dataset, year):
    """Path of one downloaded standard-level file."""
    return Path(wod_std_dir) / dataset.lower() / f'{dataset}S{year}.gz'


def index_file(wod_std_dir, dataset, year):
    """Path of one file's index (made by :func:`build`)."""
    return Path(wod_std_dir) / 'index' / f'{dataset}S{year}.npz'


# --------------------------------------------------------------------------------------------------
# download and index
# --------------------------------------------------------------------------------------------------
def download(wod_std_dir, datasets=DATASETS, years=None):
    """Download the standard-level files not already present; returns the number downloaded."""
    years = years or range(FIRST_YEAR, datetime.date.today().year + 1)
    n = 0
    for ds in datasets:
        for y in years:
            out = std_file(wod_std_dir, ds, y)
            if out.exists() and out.stat().st_size > 0:
                continue
            out.parent.mkdir(parents=True, exist_ok=True)
            url = f'{BASE}/{ds}/STD/{ds}S{y}.gz'
            tmp = out.with_suffix('.gz.part')
            try:
                urllib.request.urlretrieve(url, tmp)
            except urllib.error.HTTPError as e:
                if e.code == 404:          # no data of this kind in this year
                    continue
                raise
            tmp.rename(out)
            n += 1
            print(f'  downloaded {out.name} ({out.stat().st_size / 1e6:.1f} MB)', flush=True)
    return n


def _read_std_file(args):
    """Read one standard-level file into flat arrays (see :func:`build`)."""
    src, dst = args
    from wodpy import wod
    cast, off, z, flags, pflags = [], [0], [], {v: [] for v in _VAR_CODES}, {v: [] for v in _VAR_CODES}
    with tempfile.TemporaryDirectory() as td:   # wodpy needs an uncompressed file
        raw = Path(td) / 'wod'
        with gzip.open(src, 'rb') as fi, open(raw, 'wb') as fo:
            shutil.copyfileobj(fi, fo, 1 << 24)
        if raw.stat().st_size == 0:
            return str(dst), 0
        with open(raw) as fid:
            while True:
                p = wod.WodProfile(fid)
                zz = np.ma.filled(np.ma.asarray(p.z(), float), np.nan)
                cast.append(p.uid()); z.append(zz); off.append(off[-1] + len(zz))
                for v, code in _VAR_CODES.items():
                    i = p.var_index(code)
                    if i is None:
                        flags[v].append(np.full(len(zz), -1, np.int8)); pflags[v].append(-1)
                    else:   # -1 where the level has no value of this variable
                        flags[v].append(np.ma.filled(p.var_level_qc(i), -1).astype(np.int8))
                        pflags[v].append(p.var_profile_qc(i) or 0)
                if p.is_last_profile_in_file(fid):
                    break
    cast = np.array(cast, np.int64); order = np.argsort(cast, kind='stable')
    lens = np.diff(off); off = np.array(off, np.int64)
    take = np.concatenate([np.arange(off[k], off[k + 1]) for k in order]) if len(cast) else np.array([], int)
    out = dict(cast=cast[order], off=np.concatenate([[0], np.cumsum(lens[order])]), z=np.concatenate(z)[take].astype(np.float32))
    for v in _VAR_CODES:
        out[f'{v}flag'] = np.concatenate(flags[v])[take]
        out[f'{v}profile_flag'] = np.array(pflags[v], np.int8)[order]
    dst.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(dst.with_suffix('.part.npz'), **out)
    dst.with_suffix('.part.npz').rename(dst)
    return str(dst), len(cast)


def build(wod_std_dir, n_workers=None, rebuild=False):
    """Read every downloaded standard-level file once into ``<wod_std_dir>/index/<DATASET>S<year>.npz``:
    the casts sorted by cast number, with their levels' depths and flags (T and S; -1 = no value),
    and the whole-cast flags. Files whose index is newer are skipped unless ``rebuild``."""
    jobs = []
    for ds in DATASETS:
        for src in sorted((Path(wod_std_dir) / ds.lower()).glob(f'{ds}S*.gz')):
            dst = index_file(wod_std_dir, ds, src.name[len(ds) + 1:-3])
            if rebuild or not dst.exists() or dst.stat().st_mtime < src.stat().st_mtime:
                jobs.append((src, dst))
    print(f'indexing {len(jobs)} standard-level files', flush=True)
    n_workers = n_workers or min(8, max(1, (os.cpu_count() or 1) - 1))
    with ProcessPoolExecutor(n_workers) as ex:
        for dst, n in ex.map(_read_std_file, jobs):
            print(f'  {Path(dst).name}: {n:,} casts', flush=True)
    return len(jobs)


def check(wod_std_dir):
    """Problems with the setup (empty list if it is usable)."""
    if not wod_std_dir:
        return ['paths.wod_std_dir is not set']
    d = Path(wod_std_dir)
    if not d.is_dir():
        return [f'paths.wod_std_dir does not exist: {d}']
    src = [p for ds in DATASETS for p in (d / ds.lower()).glob(f'{ds}S*.gz')]
    idx = list((d / 'index').glob('*.npz'))
    probs = []
    if not src:
        probs.append(f'no standard-level files in {d}/<dataset>/; run: ecco-insitu wod-flags --download')
    missing = [p.name for p in src if not index_file(d, p.name.split('S')[0], p.name[len(p.name.split('S')[0]) + 1:-3]).exists()]
    if src and missing:
        probs.append(f'{len(missing)} of {len(src)} files not indexed (e.g. {missing[0]}); run: ecco-insitu wod-flags --build')
    if not idx and not probs:
        probs.append('no index; run: ecco-insitu wod-flags --build')
    return probs


# --------------------------------------------------------------------------------------------------
# applying the flags to a profile file
# --------------------------------------------------------------------------------------------------
def cast_numbers(descr):
    """WOD instrument code and cast number per profile from ``prof_descr`` (code '' and -1 if none)."""
    codes, casts = [], []
    for s in descr:
        s = s.decode(errors='replace') if isinstance(s, bytes) else str(s)
        m = _CAST.match(s.strip())
        codes.append(m.group(1) if m else ''); casts.append(int(m.group(2)) if m else -1)
    return np.array(codes), np.array(casts, np.int64)


#: Default for ``chain.wod_flags_extend_levels`` (see the module description; pending approval).
EXTEND_LEVELS = 2


def level_flags(z_std, f_std, z, extend_levels=EXTEND_LEVELS):
    """Flags on depths ``z`` from a cast's standard levels ``z_std`` (ascending) and their flags ``f_std``
    (-1 = no value; those levels are skipped): the flag at the same depth, else the larger of the two
    levels around it. Outside the cast's levels: the nearest end level's flag if the ``extend_levels``
    levels at that end are all flagged (0 = never), else 0."""
    ok = f_std >= 0
    zs, fs = z_std[ok], f_std[ok].astype(np.int16)
    out = np.zeros(len(z), np.int16)
    if len(zs) == 0:
        return out
    j = np.searchsorted(zs, z, side='left')                       # zs[j-1] < z <= zs[j]
    jc = np.minimum(j, len(zs) - 1)
    exact = np.isclose(zs[jc], z)
    inside = (j >= 1) & (j < len(zs))
    lo, hi = fs[np.clip(j - 1, 0, len(zs) - 1)], fs[jc]
    out = np.where(exact, fs[jc], np.where(inside, np.maximum(lo, hi), 0))
    n = int(extend_levels)
    if n > 0 and len(fs) >= n:
        if (fs[-n:] > 0).all():
            out = np.where(z > zs[-1] + 1e-6, fs[-1], out)
        if (fs[:n] > 0).all():
            out = np.where(z < zs[0] - 1e-6, fs[0], out)
    return out.astype(np.int16)


def _load_index(wod_std_dir, dataset, year, _cache={}):
    key = (str(wod_std_dir), dataset, year)
    if key not in _cache:
        p = index_file(wod_std_dir, dataset, year)
        _cache.clear()   # one input file (one year) at a time per worker
        _cache[key] = dict(np.load(p)) if p.exists() else None
    return _cache[key]


def flags_for(descr, yyyymmdd, depth, wod_std_dir, variables=('prof_T', 'prof_S'), extend_levels=EXTEND_LEVELS):
    """WOD flags (iPROF x iDEPTH, int16) for each variable, and counts for the log.

    Each profile is looked up in the datasets for its instrument code first, then in the others, in its
    year and the years either side (a date can move across New Year when the chain normalizes it)."""
    codes, casts = cast_numbers(descr)
    years = (np.asarray(yyyymmdd, float) // 10000).astype(int)
    depth = np.asarray(depth, float)
    out = {v: np.zeros((len(casts), len(depth)), np.int16) for v in variables}
    found = np.zeros(len(casts), bool)
    groups = {}
    for k in np.nonzero(casts > 0)[0]:
        groups.setdefault((codes[k], years[k]), []).append(k)
    for (code, y), ks in groups.items():
        ks = np.array(ks)
        first = CODE_TO_DATASETS.get(code, ())
        for ds in first + tuple(d for d in DATASETS if d not in first):
            for yy in (y, y - 1, y + 1):
                todo = ks[~found[ks]]
                if len(todo) == 0:
                    break
                ix = None
                if 1900 < yy < 2100:
                    ix = _load_index(wod_std_dir, ds, yy)
                if ix is None:
                    continue
                pos = np.searchsorted(ix['cast'], casts[todo])
                pos = np.minimum(pos, len(ix['cast']) - 1)
                hit = ix['cast'][pos] == casts[todo]
                for k, i in zip(todo[hit], pos[hit]):
                    s = slice(ix['off'][i], ix['off'][i + 1])
                    zs = ix['z'][s].astype(float)
                    for v in variables:
                        if f'{v}flag' in ix:
                            out[v][k] = level_flags(zs, ix[f'{v}flag'][s], depth, extend_levels)
                found[todo[hit]] = True
    counts = {'with_cast_number': int((casts > 0).sum()), 'found': int(found.sum())}
    return out, counts


def apply(MITprof_ds, wod_std_dir, variables=('prof_T', 'prof_S'), extend_levels=EXTEND_LEVELS):
    """Add ``prof_Vflag`` (WOD's per-level flags on the file's depths) for each variable in the file.

    Profiles from outside WOD (no ``WOD-...//<cast>`` in ``prof_descr``) are left without flags (0). If no
    profile has a WOD cast number, the dataset is returned unchanged. Prints a summary for the log.
    """
    import xarray as xr
    if 'prof_descr' not in MITprof_ds or 'prof_depth' not in MITprof_ds:
        return MITprof_ds
    vs = [v for v in variables if v in MITprof_ds]
    descr = MITprof_ds['prof_descr'].values
    if not any(_CAST.match((s.decode(errors='replace') if isinstance(s, bytes) else str(s)).strip()) for s in descr[:1000]):
        print('  WOD flags: no WOD cast numbers in prof_descr; none added')
        return MITprof_ds
    flags, c = flags_for(descr, MITprof_ds['prof_YYYYMMDD'].values, MITprof_ds['prof_depth'].values, wod_std_dir, vs, extend_levels)
    print(f"  WOD flags: {c['found']:,} of {c['with_cast_number']:,} WOD casts found in the WOD standard-level index")
    for v in vs:
        dims = MITprof_ds[v].dims
        has = MITprof_ds[v].notnull().values
        n_flag = int(((flags[v] > 0) & has).sum())
        print(f"  WOD flags: {v}: {n_flag:,} of {int(has.sum()):,} values flagged ({100 * n_flag / max(int(has.sum()), 1):.2f}%)")
        if f'{v}flag' in MITprof_ds:   # keep any flag the input already has
            flags[v] = np.maximum(flags[v], np.nan_to_num(MITprof_ds[f'{v}flag'].values, nan=0).astype(np.int16))
        MITprof_ds[f'{v}flag'] = xr.DataArray(flags[v], dims=dims, attrs=flag_attrs(v, extend_levels))
    return MITprof_ds


def flag_attrs(v, extend_levels=EXTEND_LEVELS):
    """Attributes of ``prof_Vflag`` (the chain's steps drop attributes, so the pipeline sets them again)."""
    beyond = (f"beyond the cast's WOD levels, the nearest level's flag if the {extend_levels} levels at that end are "
              f"all flagged, else 0" if extend_levels else "beyond the cast's WOD levels, 0")
    return {'long_name': f'WOD quality flag of {v[5:]} (0 = accepted)',
            'comment': ("WOD standard-level per-level flag (WOD/YEARLY/<dataset>/STD), by WOD cast number; "
                        "on depths between WOD standard levels, the larger of the two levels' flags; "
                        f"{beyond}. Step 06 zeroes the weight where > 0."),
            'wod_flags_extend_levels': int(extend_levels)}


# --------------------------------------------------------------------------------------------------
def cli(argv=None):
    """Command-line entry point (``ecco-insitu wod-flags``); see the module description.

    Returns
    -------
    int
        0 if the setup is usable (after downloading or building, if asked), else 1.
    """
    import ecco_insitu.config as config_loader
    ap = argparse.ArgumentParser(description=__doc__.split('\n\n')[0],
                                 formatter_class=argparse.RawDescriptionHelpFormatter, epilog=__doc__)
    ap.add_argument('--download', action='store_true', help='download the standard-level files not yet present')
    ap.add_argument('--build', action='store_true', help='index the downloaded files (needs wodpy)')
    ap.add_argument('--rebuild', action='store_true', help='with --build: index every file again')
    ap.add_argument('--datasets', nargs='+', default=list(DATASETS), help=f'for --download (default: {" ".join(DATASETS)})')
    ap.add_argument('--years', nargs=2, type=int, metavar=('FIRST', 'LAST'), help=f'for --download (default: {FIRST_YEAR} to this year)')
    ap.add_argument('-n', type=int, help='parallel workers for --build')
    config_loader.add_arguments(ap)
    a = ap.parse_args(argv)
    try:
        cfg = config_loader.from_args(a)
    except config_loader.ConfigError as e:
        sys.exit(f'config error: {e}')
    d = cfg['paths']['wod_std_dir']
    if not d:
        sys.exit('set paths.wod_std_dir (where the WOD standard-level files go) in config/local.yaml first')
    if a.download:
        years = range(a.years[0], a.years[1] + 1) if a.years else None
        print(f'{download(d, a.datasets, years)} files downloaded into {d}')
    if a.build or a.rebuild:
        build(d, a.n, rebuild=a.rebuild)
    probs = check(d)
    print(f'paths.wod_std_dir = {d}\nchain.apply_wod_flags = {cfg["chain"]["apply_wod_flags"]}')
    print('WOD flags: ' + ('OK' if not probs else 'NOT READY\n  ' + '\n  '.join(probs)))
    return 0 if not probs else 1


if __name__ == '__main__':
    sys.exit(cli())
