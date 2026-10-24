"""
Shared helpers for the diagnostics reports (profile counts, costs): default paths, run selection,
source/subtype parsing, geodesic bin lookup and the lat/lon bin raster used
for maps.
"""

import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree

import ecco_insitu.tools as tools  # noqa: E402  (repo-level tools.py: LLC90 grid loader, sph2cart)
import ecco_insitu.config as config_loader  # noqa: E402

# ---------------------------------------------------------------------------
# Default locations, from the chain settings (ecco_insitu/defaults.yaml + config/local.yaml):
# the same grid/bin files as ecco_insitu/pipeline.py step 02. A setting that isn't set gives None here;
# the scripts then need the matching command-line option.
# ---------------------------------------------------------------------------

_CFG = config_loader.load()


def _p(v):
    return Path(v) if v else None


RAW_ROOT         = _p(_CFG['run']['input_dir'])
PROCESSED_PARENT = _p(Path(_CFG['run']['output_base']).parent if _CFG['run']['output_base'] else None)
REPORT_PARENT    = _p(_CFG['diagnostics']['report_dir'])
GRID_DIR         = _p(_CFG['paths']['grid_dir'])
SPHERE_BIN_DIR   = _p(_CFG['paths']['sphere_bin_dir'])

# grid key -> (bin-id variable in processed files, number of bins, LLC90 id file)
GRIDS = {
    'a': ('prof_bin_id_a', 10242, 'llc090_sphere_point_n_10242_ids.bin'),
    'b': ('prof_bin_id_b',  2562, 'llc090_sphere_point_n_02562_ids.bin'),
}
VARS = ('T', 'S')

_RUN_RE = re.compile(r'^profile_files_NCEI_processed_(\d{8}_\d{6})$')
_STEP_RE = re.compile(r'__ncei_step_(\d+)$')

# ---------------------------------------------------------------------------
# Run selection
# ---------------------------------------------------------------------------

def run_root(path):
    """The timestamped run folder, given it or its yearly/ subfolder."""
    p = Path(path)
    return p.parent if p.name == 'yearly' else p


def yearly_dir(path):
    """Where a run's per-year files are: the run folder itself, or a
    <run>/yearly subfolder if the output layout has one."""
    root = run_root(path)
    return root / 'yearly' if (root / 'yearly').is_dir() else root


def run_timestamp(run_dir):
    """The ``YYYYMMDD_HHMMSS`` stamp of a run folder ``profile_files_NCEI_processed_<stamp>``
    (given it or its ``yearly/`` subfolder); the folder name itself if it doesn't follow that pattern.
    """
    name = run_root(run_dir).name
    m = _RUN_RE.match(name)
    return m.group(1) if m else name


def latest_processed_run(parent=PROCESSED_PARENT, min_step=11):
    """Newest profile_files_NCEI_processed_<ts> dir whose files all reached
    step >= min_step (i.e. a finished run with spatial scaling). Runs still in
    progress, backups and other suffixed copies are skipped.  Returns the
    folder holding the per-year files (see yearly_dir)."""
    runs = sorted((d for d in Path(parent).iterdir()
                   if d.is_dir() and _RUN_RE.match(d.name)),
                  key=lambda d: d.name, reverse=True)
    for d in runs:
        y = yearly_dir(d)
        steps = [_STEP_RE.search(f.stem) for f in y.rglob('*.nc')]
        if steps and all(m and int(m.group(1)) >= min_step for m in steps):
            return y
        print(f"  [run select] skipping {d.name}: not every file at step >= {min_step}")
    raise FileNotFoundError(f"No finished step-{min_step} run under {parent}")

# ---------------------------------------------------------------------------
# Source / subtype
# ---------------------------------------------------------------------------

# Label used when a file carries no finer qualifier
_DEFAULT_SUBTYPE = {'XBT_WOD': 'flag'}


def classify(path, root):
    """Return (source, subtype) for a profile file.

    - **source**: the top-level directory under ``root`` (CTD_WOD, PFL, ...).
    - **subtype**: the nested subdirectory if any (ITP/L2); else the filename tokens after the year
      and instrument code (PFL ``D``, XBT ``noflag``, CTD ``OSD``); else a per-source default.
    """
    rel = Path(path).relative_to(root)
    source = rel.parts[0]
    if len(rel.parts) > 2:
        return source, '/'.join(rel.parts[1:-1])
    stem = _STEP_RE.sub('', Path(path).stem)
    tokens = stem.split('_')
    yi = next((i for i, t in enumerate(tokens) if re.fullmatch(r'(19|20)\d\d', t)), None)
    rest = tokens[yi + 2:] if yi is not None else []
    return source, ('_'.join(rest) or _DEFAULT_SUBTYPE.get(source, 'all'))


def file_year(path):
    """The year in a file name (``..._1992_...``, ignoring the ``__ncei_step_N`` suffix), or -1."""
    m = re.search(r'_((?:19|20)\d\d)(?:_|$)', _STEP_RE.sub('', Path(path).stem))
    return int(m.group(1)) if m else -1

# ---------------------------------------------------------------------------
# Geodesic bins (same nearest-LLC90-cell lookup as step02)
# ---------------------------------------------------------------------------

class BinLookup:
    """lon/lat -> geodesic bin ids via the nearest LLC90 cell, as in step02."""

    def __init__(self, grid_dir=GRID_DIR, sphere_bin_dir=SPHERE_BIN_DIR):
        lon90, lat90, bathy90, X, Y, Z = tools.load_llc90_grid_step2(str(grid_dir))
        self._tree = cKDTree(np.column_stack((X.ravel(), Y.ravel(), Z.ravel())))
        self._land = (bathy90.ravel() == 0)
        self._ids = {}
        for g, (_, _, fname) in GRIDS.items():
            with open(Path(sphere_bin_dir) / fname, 'rb') as fid:
                ids = np.fromfile(fid, dtype='>f4').reshape(
                    (90, 13 * 90, 1), order='F').ravel().astype(np.int32)
            # the 2562 file marks some cells -9999; treat as unassigned (0)
            self._ids[g] = np.where(ids >= 1, ids, 0)

    def _cells(self, lon, lat):
        lon = np.asarray(lon, dtype=np.float64)
        lat = np.asarray(lat, dtype=np.float64)
        ok = np.isfinite(lon) & np.isfinite(lat) & (np.abs(lat) <= 90)
        cells = np.full(lon.shape, -1, dtype=np.int64)
        if ok.any():
            d2r = np.pi / 180.0
            xyz = np.column_stack(tools.sph2cart(lon[ok] * d2r, lat[ok] * d2r, 1.0))
            cells[ok] = self._tree.query(xyz)[1]
        return cells

    def bins(self, lon, lat):
        """{grid: int32 bin ids}; 0 where lon/lat is invalid."""
        cells = self._cells(lon, lat)
        ok = cells >= 0
        out = {}
        for g, ids in self._ids.items():
            b = np.zeros(cells.shape, dtype=np.int32)
            b[ok] = ids[cells[ok]]
            out[g] = b
        return out

    def projected_raster(self, view='global', npix=900):
        """Raster in a map projection: {grid: bin id per pixel (0 = none)},
        land mask, outside-globe mask, imshow extent.

        view: 'global' (Robinson), 'NH' / 'SH' (polar stereographic, to 50°).
        Pixels are inverse-projected to lon/lat and binned with the step 02
        lookup, so drawing a map is just per_bin[raster].
        """
        import cartopy.crs as ccrs
        pc = ccrs.PlateCarree()
        if view == 'global':
            proj = ccrs.Robinson()
            x0, x1 = proj.x_limits
            y0, y1 = proj.y_limits
            nx, ny = npix, int(npix * (y1 - y0) / (x1 - x0))
            lat_min = None
        else:
            proj = ccrs.NorthPolarStereo() if view == 'NH' else ccrs.SouthPolarStereo()
            lat_min = 50.0
            r = abs(proj.transform_point(0, lat_min if view == 'NH' else -lat_min, pc)[1])
            x0, x1, y0, y1 = -r, r, -r, r
            nx = ny = int(npix * 0.6)
        xs = np.linspace(x0, x1, nx)
        ys = np.linspace(y0, y1, ny)
        X, Y = np.meshgrid(xs, ys)
        ll = pc.transform_points(proj, X.ravel(), Y.ravel())
        lon, lat = ll[:, 0], ll[:, 1]
        outside = ~(np.isfinite(lon) & np.isfinite(lat))
        if lat_min is not None:
            outside |= (np.abs(lat) < lat_min)
        cells = self._cells(np.where(outside, 0, lon), np.where(outside, 0, lat))
        shape = X.shape
        bins = {g: np.where(outside, 0, ids[cells]).reshape(shape) for g, ids in self._ids.items()}
        land = (self._land[cells] & ~outside).reshape(shape)
        return bins, land, outside.reshape(shape), (x0, x1, y0, y1)

    def raster(self, res=0.5):
        """Lat/lon raster for maps: {grid: bin id per pixel}, land mask, extent."""
        lons = np.arange(-180 + res / 2, 180, res)
        lats = np.arange(-90 + res / 2, 90, res)
        lon2, lat2 = np.meshgrid(lons, lats)
        cells = self._cells(lon2.ravel(), lat2.ravel())
        shape = lon2.shape
        bins = {g: ids[cells].reshape(shape) for g, ids in self._ids.items()}
        land = self._land[cells].reshape(shape)
        return bins, land, (-180, 180, -90, 90)

# ---------------------------------------------------------------------------
# Index I/O
# ---------------------------------------------------------------------------

def index_paths(report_dir):
    """Paths of a profile-counts report's index files:
    ``(<report_dir>/index/raw.parquet, processed.parquet, meta.json)``.
    """
    d = Path(report_dir) / 'index'
    return d / 'raw.parquet', d / 'processed.parquet', d / 'meta.json'


def load_index(report_dir):
    """Read a report's per-profile index; returns ``(raw, processed)`` DataFrames."""
    raw_p, proc_p, _ = index_paths(report_dir)
    return pd.read_parquet(raw_p), pd.read_parquet(proc_p)


def sources_in(*dfs):
    """Sorted list of the ``source`` values found in any of the given DataFrames."""
    return sorted(set().union(*(set(df['source'].unique()) for df in dfs)))
