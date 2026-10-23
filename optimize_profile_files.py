#!/usr/bin/env python3
"""
Arrange profile files for the MITgcm profiles package, minimizing the per-tile
array size NFILESPROFMAX x NOBSGLOB, for any LLC grid and tile size.

Plain-language guide (what this does and why): docs/profile_file_optimization.md

The model reads every profile file and sends each profile to the tile it lies
in; per tile it allocates arrays of NFILESPROFMAX (files) x NOBSGLOB (the most
profiles any one file puts in any one tile) x NLEVELMAX.  Profiles are a fixed
number, so minimizing NFILESPROFMAX x NOBSGLOB is minimizing empty space.

Rules for the files:
  * one provenance per file (data.profiles weights are per file),
  * chronological, each file one continuous time span,
  * year rule (--year_rule):
      whole_years (default) files are labeled in whole calendar years (1 Jan to
                            31 Dec of the same or a later year): a file holds
                            several whole years, or, if one year alone has more
                            than NOBSGLOB profiles in some tile, one part of that
                            year (parts share the year's label, e.g.
                            XBT_WOD_001_19920101_19921231, XBT_WOD_002_19920101_19921231);
                            a file never mixes part of a year with other years
      none                  files may start and end anywhere in time
      one_year              no file spans two calendar years
  * no file puts more than NOBSGLOB profiles in any model tile.
For a given NOBSGLOB the split is greedy in time, which gives the fewest files.

Tiles: the LLC grid's 13 N x N faces, cut into t x t tiles (the per-rank array
size, sNx = sNy) and numbered as MITgcm exch2 does (facet by facet, x fastest);
13 (N/t)^2 tiles in total (117 for LLC90 with 30x30).  The only grid inputs are
N, t and the grid's cell coordinates.  Which tiles a model run blanks, and its
bathymetry, do not matter here: the model ignores profiles on land or in blank
tiles itself.

Profile positions (--position): 'interp' uses the cell the chain's step 01 put
each profile on (prof_interp_lon / lat — right when the chain ran on this same
grid); 'obs' uses prof_lon / lat and the nearest cell of this grid.  'auto'
(default) = interp for LLC90 (the chain's grid), obs otherwise.

Choosing NOBSGLOB (--nobsglob):
  an integer, or
  auto  sweep every value (steps of --step) and take the fewest files whose
        NFILESPROFMAX x NOBSGLOB is within --tolerance of the smallest found
        (--tolerance 0 = the minimum itself); --max_files caps the file count.

Usage:
    python optimize_profile_files.py <profiles_dir>                      # report only
    python optimize_profile_files.py <profiles_dir> --write <out_dir>    # and write the files
    python optimize_profile_files.py <dir> --llc 270 --grid_dir <llc270 grid>
    python optimize_profile_files.py <dir> --nobsglob 35000 --year_rule none

Outputs (with --write OUT_DIR): OUT_DIR holds only the profile files.  The
summary files go to --docs_dir (default: a 'documentation' folder beside OUT_DIR,
i.e. <run>/documentation for <run>/ecco_arrays): ecco_arrays_manifest.json /
.csv, optimize_sweep.csv, and PROFILES_SIZE_and_data.profiles.txt (NOBSGLOB /
NFILESPROFMAX to set and the profilesfiles(...) list).  Each file records where
its profiles came from (prof_yearly_file -> global attribute yearly_files,
prof_yearly_index).
"""

import argparse
import concurrent.futures
import csv
import datetime
import json
import os
from collections import defaultdict
from pathlib import Path

import netCDF4
import numpy as np
from scipy.spatial import cKDTree

# ---------------------------------------------------------------------------
# Defaults (ECCO v4r6)
# ---------------------------------------------------------------------------

DEFAULTS = {
    'llc': 90,
    'tile': 30,
    'grid_dir': None,                         # required: run() and --grid_dir (default: paths.grid_dir in the settings)
    'split_by_subdir': ['ITP'],               # these sources split by sub-folder (ITP_L2, ITP_L3)
    'year_rule': 'whole_years',
    'nobsglob': 'auto',
    'tolerance': 0.05,
    'step': 250,
}
COMPLEVEL = 4
_YEAR = 10_000_000_000                        # time key = YYYYMMDDhhmmss


# ---------------------------------------------------------------------------
# Grid and tiles
# ---------------------------------------------------------------------------

def read_llc_field(path, n=None):
    """A compact (n, 13n) big-endian float32 LLC field (e.g. XC.data)."""
    a = np.fromfile(path, dtype='>f4')
    n = n or int(round(np.sqrt(a.size / 13)))
    return a[: 13 * n * n].reshape((n, 13 * n), order='F').astype(np.float64)


def exch2_tile_numbers(n, t):
    """1-based exch2 tile number of every cell of the compact (n, 13n) grid."""
    if n % t:
        raise SystemExit(f'tile size {t} does not divide LLC size {n}')
    tile = np.zeros((n, 13 * n), int)
    j0, first = 0, 1
    for nx, ny in [(n, 3 * n), (n, 3 * n), (n, n), (3 * n, n), (3 * n, n)]:
        cols = nx * ny // n
        i, j = np.meshgrid(np.arange(nx), np.arange(ny), indexing='ij')
        f = first + (j // t) * (nx // t) + (i // t)
        tile[:, j0:j0 + cols] = f.reshape(-1, order='F').reshape((n, cols), order='F')
        j0 += cols
        first += (nx // t) * (ny // t)
    return tile


class Grid:
    def __init__(self, grid_dir, llc, tile):
        g = Path(grid_dir)
        xc = next(p for p in (g / 'no_blank' / 'XC.data', g / 'XC.data') if p.exists())
        yc = xc.with_name('YC.data')
        self.n, self.t = llc, tile
        self.lon, self.lat = read_llc_field(xc, llc), read_llc_field(yc, llc)
        self.tile = exch2_tile_numbers(llc, tile)
        self.n_tiles = int(self.tile.max())
        d2r = np.pi / 180
        lo, la = self.lon.ravel(order='F') * d2r, self.lat.ravel(order='F') * d2r
        self._tree = cKDTree(np.column_stack((np.cos(la) * np.cos(lo), np.cos(la) * np.sin(lo), np.sin(la))))
        self._tile_flat = self.tile.ravel(order='F')

    def tiles_of(self, lon, lat):
        d2r = np.pi / 180
        lo, la = np.asarray(lon, float) * d2r, np.asarray(lat, float) * d2r
        _, cell = self._tree.query(np.column_stack((np.cos(la) * np.cos(lo), np.cos(la) * np.sin(lo), np.sin(la))))
        return self._tile_flat[cell]



# ---------------------------------------------------------------------------
# Index: one row per profile
# ---------------------------------------------------------------------------

def provenance(path, root, split_by_subdir):
    rel = Path(path).relative_to(root)
    src = rel.parts[0]
    return f'{src}_{rel.parts[1]}' if src in split_by_subdir and len(rel.parts) > 2 else src


def build_index(profiles_dir, grid, position, split_by_subdir):
    files = sorted(Path(profiles_dir).rglob('*.nc'))
    rows = defaultdict(list)
    for fi, f in enumerate(files):
        with netCDF4.Dataset(f) as nc:
            nc.set_auto_mask(False)
            v = nc.variables
            n = nc.dimensions['iPROF'].size
            if n == 0:
                continue
            lon_v, lat_v = ('prof_interp_lon', 'prof_interp_lat') if position == 'interp' else ('prof_lon', 'prof_lat')
            tl = grid.tiles_of(v[lon_v][:].ravel(), v[lat_v][:].ravel())
            t = v['prof_YYYYMMDD'][:].astype(np.int64) * 1_000_000 + v['prof_HHMMSS'][:].astype(np.int64)
        rows[provenance(f, profiles_dir, split_by_subdir)].append((np.full(n, fi), np.arange(n), t, tl))
    index = {}
    for p, parts in rows.items():
        fi, li, t, tl = (np.concatenate(x) for x in zip(*parts))
        o = np.lexsort((li, fi, t))
        index[p] = {'file': fi[o], 'local': li[o], 'time': t[o], 'tile': tl[o]}
    return files, index


# ---------------------------------------------------------------------------
# Splitting
# ---------------------------------------------------------------------------

def _occ(tile):
    return [np.flatnonzero(tile == k) for k in np.unique(tile)]


def split_cap(tile, cap, occ=None):
    """Time-contiguous chunks with <= cap profiles per tile (greedy)."""
    occ = [o for o in (occ if occ is not None else _occ(tile)) if len(o) > cap]
    n, s, out = len(tile), 0, []
    while s < n:
        cut = n
        for o in occ:
            i = np.searchsorted(o, s)
            if i + cap < len(o):
                cut = min(cut, int(o[i + cap]))
        out.append((s, cut))
        s = cut
    return out


class Splitter:
    """Precomputes per-year counts so a sweep over caps is fast."""

    def __init__(self, d, n_tiles, year_rule):
        self.tile, self.rule = d['tile'], year_rule
        year = d['time'] // _YEAR
        b = np.flatnonzero(np.diff(year)) + 1
        self.years = [(int(year[a]), int(a), int(z)) for a, z in
                      zip(np.concatenate(([0], b)), np.concatenate((b, [len(year)])))]
        self.year_counts = [np.bincount(self.tile[a:z], minlength=n_tiles + 1)[1:] for _, a, z in self.years]
        self.min_cap = 1
        self._occ_all = _occ(self.tile) if year_rule == 'none' else None
        self._occ_year = ([_occ(self.tile[a:z]) for _, a, z in self.years]
                          if year_rule in ('one_year', 'whole_years') else None)

    def chunks(self, cap):
        if self.rule == 'none':
            return split_cap(self.tile, cap, self._occ_all)
        if self.rule == 'one_year':
            out = []
            for (y, a, z), occ in zip(self.years, self._occ_year):
                out += [(a + x, a + w) for x, w in split_cap(self.tile[a:z], cap, occ)]
            return out
        # whole_years: as many whole years per file as fit under the cap; a year
        # that alone exceeds the cap gets files of its own (parts of that year)
        out, cur, start = [], None, None
        for (y, a, z), c, occ in zip(self.years, self.year_counts, self._occ_year):
            if c.max() > cap:
                if cur is not None:
                    out.append((start, a)); cur = None
                out += [(a + x, a + w) for x, w in split_cap(self.tile[a:z], cap, occ)]
                continue
            if cur is not None and (cur + c).max() > cap:
                out.append((start, a)); cur = None
            if cur is None:
                cur, start = c.copy(), a
            else:
                cur += c
        if cur is not None:
            out.append((start, self.years[-1][2]))
        return out


def sweep(splitters, step, n_used, n_tiles):
    lo = max(step, max(s.min_cap for s in splitters.values()))
    hi = max(max(int(c.max()) for c in s.year_counts) if s.rule == 'one_year' else
             int(max(np.sum(s.year_counts, axis=0).max(), 1)) for s in splitters.values())
    caps = sorted(set([lo] + list(range(((lo + step - 1) // step) * step, hi + step, step)) + [hi]))
    rows = []
    for cap in caps:
        F = sum(len(s.chunks(cap)) for s in splitters.values())
        P = cap * F
        rows.append({'nobsglob': cap, 'files': F, 'product': P,
                     'empty_pct': round(100 * (1 - n_used / (n_tiles * P)), 2)})
    return rows


def choose(rows, nobsglob, tolerance, max_files):
    if nobsglob != 'auto':
        cap = int(nobsglob)
        return next((r for r in rows if r['nobsglob'] == cap), None) or {'nobsglob': cap}
    ok = [r for r in rows if max_files is None or r['files'] <= max_files]
    if not ok:
        raise SystemExit(f'no NOBSGLOB gives <= {max_files} files')
    best = min(r['product'] for r in ok)
    near = [r for r in ok if r['product'] <= best * (1 + tolerance)]
    return min(near, key=lambda r: (r['files'], r['product']))


# ---------------------------------------------------------------------------
# Writing
# ---------------------------------------------------------------------------

def _fill_for(name, var):
    if name.endswith('weight'):
        return 0
    if var.dtype.kind == 'S':
        return b''
    fv = getattr(var, '_FillValue', None)
    if fv is not None:
        return fv
    return -9999 if var.dtype.kind in 'fi' else 0


def _schema(files, file_ids):
    schema, dims, gattrs = {}, {}, None
    for fi in file_ids:
        with netCDF4.Dataset(files[fi]) as nc:
            if gattrs is None:
                gattrs = {k: nc.getncattr(k) for k in nc.ncattrs()}
            for d, dim in nc.dimensions.items():
                if d != 'iPROF':
                    dims[d] = max(dims.get(d, 0), len(dim))
            for name, var in nc.variables.items():
                if name not in schema:
                    schema[name] = {'dims': var.dimensions, 'dtype': var.dtype,
                                    'attrs': {k: var.getncattr(k) for k in var.ncattrs()},
                                    'fill': _fill_for(name, var)}
    return schema, dims, gattrs or {}


def write_file(args):
    path, files, prov, k, rows, extra_attrs = args
    fi, li = rows['file'], rows['local']
    n = len(fi)
    file_ids = list(dict.fromkeys(fi.tolist()))
    schema, dims, gattrs = _schema(files, file_ids)
    tmp = path.with_suffix('.tmp.nc')
    with netCDF4.Dataset(tmp, 'w') as out:
        gattrs.update({'title': f'ECCO-ready profile file: {prov}, part {k}', 'provenance': prov,
                       'n_profiles': n, 'yearly_files': json.dumps([files[i].name for i in file_ids]),
                       'created_by': 'optimize_profile_files.py',
                       'created': datetime.datetime.now().isoformat(timespec='seconds'), **extra_attrs})
        out.setncatts({a: v for a, v in gattrs.items() if a != '_NCProperties'})
        out.createDimension('iPROF', n)
        for d, L in dims.items():
            out.createDimension(d, L)
        local_file = {f: j for j, f in enumerate(file_ids)}
        for name, sc in schema.items():
            fill_kw = {'fill_value': sc['attrs']['_FillValue']} if '_FillValue' in sc['attrs'] else {}
            var = out.createVariable(name, sc['dtype'], sc['dims'], zlib=True, complevel=COMPLEVEL, **fill_kw)
            var.setncatts({a: v for a, v in sc['attrs'].items() if a != '_FillValue'})
            if 'iPROF' not in sc['dims']:
                with netCDF4.Dataset(files[file_ids[0]]) as nc:
                    var[:] = nc.variables[name][:]
                continue
            ax = sc['dims'].index('iPROF')
            shape = tuple(n if d == 'iPROF' else len(out.dimensions[d]) for d in sc['dims'])
            data = np.full(shape, sc['fill'], dtype=sc['dtype'])
            for f in file_ids:
                sel = np.flatnonzero(fi == f)
                with netCDF4.Dataset(files[f]) as nc:
                    if name not in nc.variables:
                        continue
                    nc.set_auto_mask(False)
                    src = np.take(nc.variables[name][:], li[sel], axis=ax)
                idx = [slice(None)] * len(shape)
                idx[ax] = sel
                for d, (Ls, Lo) in enumerate(zip(src.shape, shape)):
                    if d != ax and Ls < Lo:
                        idx[d] = slice(0, Ls)
                data[tuple(idx)] = src
            var[:] = data
        var = out.createVariable('prof_yearly_file', 'i2', ('iPROF',), zlib=True)
        var.long_name = 'index into global attribute yearly_files'
        var[:] = np.array([local_file[f] for f in fi], np.int16)
        var = out.createVariable('prof_yearly_index', 'i4', ('iPROF',), zlib=True)
        var.long_name = 'profile index (iPROF) within that yearly file'
        var[:] = li.astype(np.int32)
    tmp.replace(path)
    return path.name


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def run(profiles_dir, out_dir=None, llc=DEFAULTS['llc'], tile=DEFAULTS['tile'], grid_dir=DEFAULTS['grid_dir'],
        position='auto', split_by_subdir=None, year_rule=DEFAULTS['year_rule'],
        nobsglob=DEFAULTS['nobsglob'], tolerance=DEFAULTS['tolerance'], max_files=None, step=DEFAULTS['step'],
        n_workers=None, docs_dir=None):
    profiles_dir = Path(profiles_dir)
    split_by_subdir = set(DEFAULTS['split_by_subdir'] if split_by_subdir is None else split_by_subdir)
    if grid_dir is None:
        raise ValueError('grid_dir is required (the folder with XC.data, YC.data)')
    grid = Grid(grid_dir, llc, tile)
    if position == 'auto':
        position = 'interp' if llc == 90 else 'obs'
    print(f'\n{"=" * 64}\nOptimize profile files: {profiles_dir}')
    print(f'  grid LLC{llc}, {tile}x{tile} tiles: {grid.n_tiles} tiles')
    print(f'  positions: {position}; year rule: {year_rule}; one provenance per file '
          f'(split by sub-folder: {sorted(split_by_subdir) or "none"})')

    files, index = build_index(profiles_dir, grid, position, split_by_subdir)
    n_total = sum(len(d['time']) for d in index.values())
    if not n_total:
        raise SystemExit(f'No profiles found under {profiles_dir}')
    n_used = n_total
    splitters = {p: Splitter(d, grid.n_tiles, year_rule) for p, d in index.items()}
    lower = sum(int(np.sum(s.year_counts, axis=0).max()) for s in splitters.values())
    print(f'  {n_total:,} profiles, {len(index)} provenances, {len(files)} input files')
    print(f'  lower bound on NFILESPROFMAX x NOBSGLOB (sum of each provenance\'s busiest tile): {lower:,}')
    if year_rule == 'whole_years':
        dens = max(int(c.max()) for s in splitters.values() for c in s.year_counts)
        print(f'  densest provenance-year in one tile: {dens:,} (below that, years are split into parts)')

    rows = sweep(splitters, step, n_used, grid.n_tiles)
    pick = choose(rows, nobsglob, tolerance, max_files)
    cap = pick['nobsglob']
    best = min(rows, key=lambda r: r['product'])
    # trade-off table: for each file count (fewest first, up to SHOW_FILES), the best NOBSGLOB
    by_F = {}
    for r in rows:
        if r['files'] not in by_F or r['product'] < by_F[r['files']]['product']:
            by_F[r['files']] = r
    SHOW_FILES = max(60, pick.get('files', 0))
    print(f'\n  trade-off (best NOBSGLOB for each number of files; full sweep in optimize_sweep.csv):')
    print(f'  {"files":>6} {"NOBSGLOB":>9} {"files x NOBSGLOB":>16} {"empty %":>8}')
    for F_ in sorted(by_F):
        if F_ > SHOW_FILES:
            break
        r = by_F[F_]
        mark = ' <- chosen' if r['nobsglob'] == pick['nobsglob'] else (' <- minimum' if r is best else '')
        print(f"  {F_:>6} {r['nobsglob']:>9,} {r['product']:>16,} {r['empty_pct']:>7.1f}%{mark}")
    print(f"\n  chosen: NOBSGLOB = {cap:,}"
          + (f", {pick['files']} files, product {pick['product']:,} ({100 * (pick['product'] / best['product'] - 1):.1f}% "
             f"above the minimum {best['product']:,})" if 'files' in pick else ''))

    plan, manifest = [], []
    for prov in sorted(index):
        d = index[prov]
        for k, (a, b) in enumerate(splitters[prov].chunks(cap), 1):
            rows_k = {key: arr[a:b] for key, arr in d.items()}
            busiest = int(np.bincount(rows_k['tile']).max())
            if year_rule == 'whole_years':
                t0, t1 = f"{rows_k['time'][0] // _YEAR}0101", f"{rows_k['time'][-1] // _YEAR}1231"
            else:
                t0, t1 = str(rows_k['time'][0])[:8], str(rows_k['time'][-1])[:8]
            name = f'{prov}_{k:03d}_{t0}_{t1}.nc'
            attrs = {'time_start': t0, 'time_end': t1, 'first_profile_time': str(rows_k['time'][0]),
                     'last_profile_time': str(rows_k['time'][-1]), 'nobsglob_cap': cap,
                     'busiest_tile_count': busiest, 'year_rule': year_rule}
            plan.append((None, files, prov, k, rows_k, attrs, name))
            manifest.append({'file': name, 'provenance': prov, 'part': k, 'time_start': t0, 'time_end': t1,
                             'n_profiles': b - a, 'busiest_tile': busiest})
    F = len(manifest)
    used = max(m['busiest_tile'] for m in manifest)
    print(f'  -> {F} files; busiest tile in any file {used:,}; set NOBSGLOB = {used:,}, NFILESPROFMAX = {F} '
          f'(product {F * used:,}, {100 * (1 - n_used / (grid.n_tiles * F * used)):.1f}% empty over {grid.n_tiles} tiles)')
    for m in manifest:
        print(f"     {m['file']:42s} {m['n_profiles']:>10,} profiles, busiest tile {m['busiest_tile']:,}")
    if out_dir is None:
        print('  (report only — pass --write <out_dir> to write the files)\n' + '=' * 64)
        return manifest

    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    docs_dir = Path(docs_dir) if docs_dir else out_dir.parent / 'documentation'
    docs_dir.mkdir(parents=True, exist_ok=True)
    with open(docs_dir / 'optimize_sweep.csv', 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
    jobs = [(out_dir / name, fl, prov, k, rk, attrs) for _, fl, prov, k, rk, attrs, name in plan]
    n_workers = n_workers or min(4, max(1, (os.cpu_count() or 2) - 1))
    with concurrent.futures.ProcessPoolExecutor(max_workers=n_workers) as ex:
        for name in ex.map(write_file, jobs):
            print(f'  wrote {name}', flush=True)
    meta = {'profiles_dir': str(profiles_dir), 'arrays_dir': str(out_dir), 'llc': llc, 'tile': tile, 'grid_dir': str(grid_dir),
            'n_tiles': grid.n_tiles, 'position': position,
            'split_by_subdir': sorted(split_by_subdir), 'year_rule': year_rule,
            'nobsglob_cap': cap, 'nobsglob_needed': used, 'nfilesprofmax_needed': F,
            'total_profiles': n_total, 'files': manifest}
    (docs_dir / 'ecco_arrays_manifest.json').write_text(json.dumps(meta, indent=2))
    with open(docs_dir / 'ecco_arrays_manifest.csv', 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=list(manifest[0])); w.writeheader(); w.writerows(manifest)
    lines = ['C PROFILES_SIZE.h — values needed for these files',
             f'      PARAMETER ( NOBSGLOB = {used} )', f'      PARAMETER ( NFILESPROFMAX = {F} )', '',
             '# data.profiles — file list (names without .nc)']
    lines += [f" profilesfiles({i})= '{m['file'][:-3]}'," for i, m in enumerate(manifest, 1)]
    (docs_dir / 'PROFILES_SIZE_and_data.profiles.txt').write_text('\n'.join(lines) + '\n')
    print(f'  done: {F} files in {out_dir}; summary files in {docs_dir}\n' + '=' * 64)
    return manifest


def main():
    ap = argparse.ArgumentParser(description='Arrange profile files to minimize NFILESPROFMAX x NOBSGLOB.',
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    ap.add_argument('profiles_dir', help='folder of input profile files (searched recursively)')
    ap.add_argument('--write', metavar='OUT_DIR', default=None, help='write the files here (default: report only)')
    ap.add_argument('--docs_dir', default=None,
                    help="summary files (manifest, sweep, PROFILES_SIZE text); default: 'documentation' beside OUT_DIR")
    ap.add_argument('--llc', type=int, default=DEFAULTS['llc'], help='LLC grid size N')
    ap.add_argument('--tile', type=int, default=DEFAULTS['tile'], help='per-rank tile size (cells per side, sNx = sNy)')
    ap.add_argument('--grid_dir', default=None,
                    help='folder with [no_blank/]XC.data, YC.data (default: paths.grid_dir in config/local.yaml)')
    ap.add_argument('--position', choices=['auto', 'interp', 'obs'], default='auto')
    ap.add_argument('--split_by_subdir', default=','.join(DEFAULTS['split_by_subdir']),
                    help='sources whose sub-folders are separate provenances (comma list, empty for none)')
    ap.add_argument('--year_rule', choices=['whole_years', 'none', 'one_year'], default=DEFAULTS['year_rule'])
    ap.add_argument('--nobsglob', default=DEFAULTS['nobsglob'], help="integer, or 'auto'")
    ap.add_argument('--tolerance', type=float, default=DEFAULTS['tolerance'],
                    help='auto: fewest files within this fraction of the minimum product (0 = minimum)')
    ap.add_argument('--max_files', type=int, default=None, help='auto: at most this many files')
    ap.add_argument('--step', type=int, default=DEFAULTS['step'], help='auto: NOBSGLOB sweep step')
    ap.add_argument('-n', '--n_workers', type=int, default=None)
    a = ap.parse_args()
    if a.grid_dir is None:
        import config_loader
        a.grid_dir = config_loader.load()['paths']['grid_dir']
        if a.grid_dir is None:
            ap.error('give --grid_dir, or set paths.grid_dir in config/local.yaml')
    run(a.profiles_dir, a.write, a.llc, a.tile, a.grid_dir, a.position,
        [s for s in a.split_by_subdir.split(',') if s], a.year_rule, a.nobsglob, a.tolerance,
        a.max_files, a.step, a.n_workers, a.docs_dir)


if __name__ == '__main__':
    main()
