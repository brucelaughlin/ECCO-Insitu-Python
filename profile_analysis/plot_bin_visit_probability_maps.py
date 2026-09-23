"""
Bin visit probability maps for all data classes.

For each data class, produces:
  - Geodesic polygon maps (native bin polygons, Robinson projection)
  - EASE-Grid 2.0 maps (nearest-neighbour interpolated, 36 km, EPSG:6933)

Both as: all-years combined + per-year (from auto-detected start year to Dec 2025).
One file per threshold (>=1, >=2, >=3 profiles/bin/month).

Output tree:
  plots/<class_label>/bin_visit_probability_maps/
  plots/<class_label>/bin_visit_probability_maps_ease2/
  plots/<class_label>/bin_visit_probability_maps_ease2/per_year/

Probability denominator:
  - All-years: total calendar months spanned by record (first to last, inclusive)
  - Per-year:  12 (months in the year); truncated at Dec 2025

Run time note: LLC90 grid and geodesic polygons are built once and reused across
all data classes. Bin assignment is re-run per class (no precomputed bin IDs).
"""

import sys
import re
import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.cm as cm
import matplotlib.collections as mcollections
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import cartopy.io.shapereader as shpreader
import shapely.ops
from pathlib import Path
from collections import defaultdict
from scipy.interpolate import griddata
from scipy.spatial import KDTree, ConvexHull
from shapely.geometry import Polygon as ShapelyPolygon
from pyproj import Transformer

sys.path.append('/Users/brucel/ecco/yip/ECCO-Insitu-Python')
import tools
from tools import sph2cart

# ==============================================================================
# Switches
# ==============================================================================

RUN_POLYGON = True   # geodesic polygon maps
RUN_EASE2   = True   # EASE-Grid 2.0 maps

# ==============================================================================
# Data class definitions
# Each entry: label (used for output folder) -> input directory Path
# ==============================================================================

INTERP_ROOT = Path('/Users/brucel/ecco/yip/profile_data/Interp_Profiles')

DATA_CLASSES = {
    'CTD_WOD':         INTERP_ROOT / 'CTD_WOD',
    'GLD_WOD':         INTERP_ROOT / 'GLD_WOD',
    'ITP_L2':          INTERP_ROOT / 'ITP/L2',
    'ITP_L3':          INTERP_ROOT / 'ITP/L3',
    'MEOP':            INTERP_ROOT / 'MEOP',
    'MRB_WOD':         INTERP_ROOT / 'MRB_WOD',
    'MRB_WOD_36depth': INTERP_ROOT / 'MRB_WOD/36_depths',
    'MRB_WOD_97depth': INTERP_ROOT / 'MRB_WOD/97_depths',
    'PFL':             INTERP_ROOT / 'PFL',
    'PFL_BGC':         INTERP_ROOT / 'PFL_BGC',
    'Samoa':           INTERP_ROOT / 'Samoa',
    'SOCAT':           INTERP_ROOT / 'SOCAT',
    'XBT_WOD':         INTERP_ROOT / 'XBT_WOD',
}

# ==============================================================================
# Paths
# ==============================================================================

sphere_bin_dir = '/Users/brucel/ecco/yip/sample_data/ecco-insitu/sweet_gdrive/grid_llc90/sphere_point_distribution'
grid_dir       = '/Users/brucel/ecco/yip/sample_data/ecco-insitu/sweet_gdrive/grid_llc90'

geodesic_file_dict = {
    'a': '/Users/brucel/ecco/yip/sample_data/ecco-insitu/sweet_gdrive/geodesic/10242_bin_locations.csv',
    'b': '/Users/brucel/ecco/yip/sample_data/ecco-insitu/sweet_gdrive/geodesic/02562_bin_locations.csv',
}
bin_labels = {'a': '10242 bins', 'b': '02562 bins'}

PLOTS_ROOT = Path('/Users/brucel/ecco/yip/profile_data/z_profile_file_analysis/plots')

# EASE-Grid 2.0 parameters
EASE2_NCOLS  = 964
EASE2_NROWS  = 406
EASE2_CELL_M = 36032.22
EASE2_X0     = -17367530.44
EASE2_Y0     =   7314540.83

# ==============================================================================
# Grid setup (built once, reused across all data classes)
# ==============================================================================

def load_bin_arrays(sphere_bin_dir):
    import os
    tile_shape = [90, 13*90, 1]
    bin_arrays = {}
    for key, fname in [('a', 'llc090_sphere_point_n_10242_ids.bin'),
                       ('b', 'llc090_sphere_point_n_02562_ids.bin')]:
        fpath = os.path.join(sphere_bin_dir, fname)
        with open(fpath, 'rb') as fid:
            arr = np.fromfile(fid, dtype='>f4')
        arr = arr.reshape((tile_shape[0], tile_shape[1], tile_shape[2]), order='F')
        bin_arrays[key] = arr.ravel()
    return bin_arrays


def build_llc90_griddata_inputs(grid_dir):
    lon_90, lat_90, bathy_90, X_90, Y_90, Z_90 = tools.load_llc90_grid_step2(grid_dir)
    xyz_grid = np.column_stack((X_90.ravel(), Y_90.ravel(), Z_90.ravel()))
    return xyz_grid, np.arange(X_90.size)


def assign_bin_ids_from_latlon(prof_lon, prof_lat, xyz_grid, grid_indices, bin_arrays):
    deg2rad = np.pi / 180.0
    valid_mask = np.asarray(
        tools.sph2cart_returnValidMaskOnly(prof_lon * deg2rad, prof_lat * deg2rad, 1),
        dtype=bool)
    x, y, z = sph2cart(prof_lon[valid_mask] * deg2rad, prof_lat[valid_mask] * deg2rad, 1)
    xyz_prof = np.column_stack((np.asarray(x), np.asarray(y), np.asarray(z)))
    llc_indices = griddata(xyz_grid, grid_indices, xyz_prof, method='nearest').astype(int)
    result = {}
    for key, bin_arr in bin_arrays.items():
        ids = np.full(len(prof_lon), np.nan)
        ids[valid_mask] = bin_arr[llc_indices]
        result[key] = ids
    return result


def load_geodesic_centres(geodesic_file):
    df = pd.read_csv(geodesic_file, header=None, names=['lon', 'lat'])
    return df['lon'].values, df['lat'].values


def build_land_polygon():
    shpfilename = shpreader.natural_earth(resolution='10m', category='physical', name='land')
    reader = shpreader.Reader(shpfilename)
    return shapely.ops.unary_union([rec.geometry for rec in reader.records()])


def build_geodesic_polygons(geodesic_file, land_polygon, angular_precision=0.5,
                            min_ocean_fraction=0.1):
    bin_lons, bin_lats = load_geodesic_centres(geodesic_file)
    xyz = np.stack(sph2cart(np.radians(bin_lons), np.radians(bin_lats), 1), axis=-1)
    tree = KDTree(xyz)

    art_lons = np.arange(-180, 180, angular_precision)
    art_lats = np.arange(-90, 90, angular_precision)
    lon_grid, lat_grid = np.meshgrid(art_lons, art_lats)
    art_xyz = np.stack(sph2cart(np.radians(lon_grid.ravel()),
                                np.radians(lat_grid.ravel()), 1), axis=-1)
    _, bin_assignments = tree.query(art_xyz)
    bin_assignments = bin_assignments.reshape(lon_grid.shape)

    polygons = {}
    for bid in range(len(bin_lons)):
        mask = bin_assignments == bid
        if not mask.any():
            continue
        pts = np.column_stack((lon_grid[mask], lat_grid[mask]))
        if len(pts) < 3:
            continue
        if pts[:, 0].max() - pts[:, 0].min() > 180:
            continue
        hull_verts = pts[ConvexHull(pts).vertices]
        bin_shapely = ShapelyPolygon(hull_verts)
        bin_ocean = bin_shapely.difference(land_polygon)
        if bin_ocean.is_empty:
            continue
        if bin_ocean.area / bin_shapely.area < min_ocean_fraction:
            continue
        if bin_ocean.geom_type == 'MultiPolygon':
            bin_ocean = max(bin_ocean.geoms, key=lambda g: g.area)
        polygons[bid] = np.array(bin_ocean.exterior.coords)
    return polygons


def build_ease2_lonlat():
    cols = np.arange(EASE2_NCOLS)
    rows = np.arange(EASE2_NROWS)
    x_centres = EASE2_X0 + (cols + 0.5) * EASE2_CELL_M
    y_centres = EASE2_Y0 - (rows + 0.5) * EASE2_CELL_M
    x_grid, y_grid = np.meshgrid(x_centres, y_centres)
    transformer = Transformer.from_crs("EPSG:6933", "EPSG:4326", always_xy=True)
    lon_flat, lat_flat = transformer.transform(x_grid.ravel(), y_grid.ravel())
    return lon_flat.reshape(x_grid.shape), lat_flat.reshape(x_grid.shape), x_grid, y_grid


def build_ease2_bin_mapping(ease2_lon, ease2_lat, geodesic_file):
    bin_lons, bin_lats = load_geodesic_centres(geodesic_file)
    bin_xyz = np.stack(sph2cart(np.radians(bin_lons), np.radians(bin_lats), 1), axis=-1)
    tree = KDTree(bin_xyz)
    ease_xyz = np.stack(sph2cart(np.radians(ease2_lon.ravel()),
                                 np.radians(ease2_lat.ravel()), 1), axis=-1)
    _, indices = tree.query(ease_xyz)
    return indices  # 0-based

# ==============================================================================
# Data collection
# ==============================================================================

def collect_records(file_list, bin_arrays, xyz_grid, grid_indices):
    records = []
    for fpath in sorted(file_list):
        try:
            ds = xr.open_dataset(fpath)
        except Exception as e:
            print(f"  skipping {fpath.name}: {e}")
            continue
        if ds.sizes.get('iPROF', 0) == 0:
            ds.close()
            continue

        yyyymmdd = ds['prof_YYYYMMDD'].values.astype(int)
        months = (yyyymmdd % 10000) // 100
        years  = yyyymmdd // 10000

        prof_lon = ds['prof_lon'].values.astype(float)
        prof_lat = ds['prof_lat'].values.astype(float)
        assigned = assign_bin_ids_from_latlon(prof_lon, prof_lat, xyz_grid, grid_indices, bin_arrays)

        for i in range(len(months)):
            records.append({
                'year':     int(years[i]),
                'month':    int(months[i]),
                'bin_id_a': float(assigned['a'][i]),
                'bin_id_b': float(assigned['b'][i]),
            })
        ds.close()
        print(f"    {fpath.name}: {len(months)} profiles")
    return records

# ==============================================================================
# Probability computation
# ==============================================================================

def compute_visit_probabilities(records, bin_key):
    ym_bin_counts = defaultdict(lambda: defaultdict(int))
    for r in records:
        bid = r[f'bin_id_{bin_key}']
        if np.isnan(bid) or bid <= 0:
            continue
        ym_bin_counts[(r['year'], r['month'])][int(bid)] += 1

    all_ym = sorted(ym_bin_counts.keys())
    if not all_ym:
        return {}
    year_min, month_min = all_ym[0]
    year_max, month_max = all_ym[-1]
    total_months = (year_max - year_min) * 12 + (month_max - month_min) + 1

    bin_month_counts = defaultdict(list)
    for ym, bin_count_dict in ym_bin_counts.items():
        for bid, cnt in bin_count_dict.items():
            bin_month_counts[bid].append(cnt)

    return {bid: {
        '1+': sum(1 for c in counts if c >= 1) / total_months,
        '2+': sum(1 for c in counts if c >= 2) / total_months,
        '3+': sum(1 for c in counts if c >= 3) / total_months,
    } for bid, counts in bin_month_counts.items()}


def compute_visit_probabilities_year(records, bin_key, year):
    ym_bin_counts = defaultdict(lambda: defaultdict(int))
    for r in records:
        if r['year'] != year or (r['year'], r['month']) > (2025, 12):
            continue
        bid = r[f'bin_id_{bin_key}']
        if np.isnan(bid) or bid <= 0:
            continue
        ym_bin_counts[(r['year'], r['month'])][int(bid)] += 1

    bin_month_counts = defaultdict(list)
    for ym, bcd in ym_bin_counts.items():
        for bid, cnt in bcd.items():
            bin_month_counts[bid].append(cnt)

    return {bid: {
        '1+': sum(1 for c in counts if c >= 1) / 12,
        '2+': sum(1 for c in counts if c >= 2) / 12,
        '3+': sum(1 for c in counts if c >= 3) / 12,
    } for bid, counts in bin_month_counts.items()}


def detect_start_year(records, bin_key, threshold_fraction=0.80):
    ym_bin_counts = defaultdict(lambda: defaultdict(int))
    for r in records:
        if (r['year'], r['month']) > (2025, 12):
            continue
        bid = r[f'bin_id_{bin_key}']
        if np.isnan(bid) or bid <= 0:
            continue
        ym_bin_counts[(r['year'], r['month'])][int(bid)] += 1

    all_years = sorted({ym[0] for ym in ym_bin_counts})
    if not all_years:
        return None

    def bins_in_year(y):
        bins = set()
        for m in range(1, 13):
            bins.update(ym_bin_counts.get((y, m), {}).keys())
        return len(bins)

    yearly_counts = {y: bins_in_year(y) for y in all_years}
    cutoff = threshold_fraction * np.median(list(yearly_counts.values()))
    for y in all_years:
        if yearly_counts[y] >= cutoff:
            return y
    return all_years[-1]

# ==============================================================================
# Plotting helpers
# ==============================================================================

CMAP = cm.RdYlGn_r
NORM = mcolors.TwoSlopeNorm(vcenter=0.5, vmin=0, vmax=1)
THRESHOLDS = ['1+', '2+', '3+']


def _thresh_label(thresh):
    return thresh.replace('+', 'plus')


def plot_polygon_maps(probs, polygons, bin_label, class_label, title, out_dir, file_stem):
    bin_ids = sorted(probs.keys())
    id_offset = 1 if (bin_ids and max(bin_ids) >= len(polygons) + 1) else 0

    for thresh in THRESHOLDS:
        verts, values = [], []
        for bid in bin_ids:
            pk = bid - id_offset
            if pk not in polygons:
                continue
            p = probs[bid][thresh]
            if p == 0:
                continue
            verts.append(polygons[pk])
            values.append(p)

        fig, ax = plt.subplots(1, 1, figsize=(14, 7),
                               subplot_kw={'projection': ccrs.Robinson()})
        ax.set_global()
        ax.add_feature(cfeature.LAND, facecolor='lightgray', zorder=1)
        ax.add_feature(cfeature.COASTLINE, linewidth=0.4, zorder=2)

        if verts:
            patches = [plt.Polygon(v) for v in verts]
            pc = mcollections.PatchCollection(patches, cmap=CMAP, norm=NORM,
                                              transform=ccrs.PlateCarree(),
                                              linewidth=0, zorder=3)
            pc.set_array(np.array(values))
            ax.add_collection(pc)
            cbar = plt.colorbar(pc, ax=ax, orientation='vertical',
                                fraction=0.02, pad=0.02)
            cbar.set_label('probability', fontsize=9)
            cbar.ax.tick_params(labelsize=8)

        ax.set_title(f'{title}\n{bin_label} (profiles assigned via LLC90),'
                     f' plotted as geodesic bin polygons'
                     f'\n>= {thresh[0]} profile(s) per bin per month', fontsize=10)
        out_dir.mkdir(parents=True, exist_ok=True)
        out = out_dir / f'{file_stem}_{_thresh_label(thresh)}.png'
        fig.tight_layout()
        fig.savefig(out, dpi=150)
        plt.close(fig)
        print(f"    saved: {out.name}")


def probs_to_ease2_grid(probs, ease2_bin_indices, n_bins):
    grids = {}
    for thresh in THRESHOLDS:
        prob_by_bin = np.full(n_bins, np.nan)
        for bid_1based, p in probs.items():
            idx = bid_1based - 1
            if 0 <= idx < n_bins:
                prob_by_bin[idx] = p[thresh]
        flat = prob_by_bin[ease2_bin_indices]
        flat[flat == 0] = np.nan
        grids[thresh] = flat.reshape(EASE2_NROWS, EASE2_NCOLS)
    return grids


def plot_ease2_maps(prob_grids, ease2_x, ease2_y, bin_label, class_label,
                   title, out_dir, file_stem):
    ease2_crs = ccrs.epsg(6933)
    for thresh in THRESHOLDS:
        fig, ax = plt.subplots(1, 1, figsize=(14, 7),
                               subplot_kw={'projection': ease2_crs})
        ax.set_global()
        ax.add_feature(cfeature.LAND, facecolor='lightgray', zorder=1)
        ax.add_feature(cfeature.COASTLINE, linewidth=0.4, zorder=2)

        masked = np.ma.masked_invalid(prob_grids[thresh])
        pm = ax.pcolormesh(ease2_x, ease2_y, masked, cmap=CMAP, norm=NORM,
                           transform=ease2_crs, zorder=3, shading='auto')

        cbar = plt.colorbar(pm, ax=ax, orientation='vertical', fraction=0.02, pad=0.02)
        cbar.set_label('probability', fontsize=9)
        cbar.ax.tick_params(labelsize=8)

        ax.set_title(f'{title}\n{bin_label} (profiles assigned via LLC90),'
                     f' interpolated to EASE-Grid 2.0 (36 km, EPSG:6933)'
                     f'\n>= {thresh[0]} profile(s) per bin per month', fontsize=10)
        out_dir.mkdir(parents=True, exist_ok=True)
        out = out_dir / f'{file_stem}_{_thresh_label(thresh)}.png'
        fig.tight_layout()
        fig.savefig(out, dpi=150)
        plt.close(fig)
        print(f"    saved: {out.name}")

# ==============================================================================
# Main
# ==============================================================================

print("Building shared grids...")
bin_arrays = load_bin_arrays(sphere_bin_dir)
xyz_grid, grid_indices = build_llc90_griddata_inputs(grid_dir)

if RUN_POLYGON:
    print("Building land polygon...")
    land_polygon = build_land_polygon()
    print("Building geodesic polygons...")
    geodesic_polygons = {}
    for bin_key in ('a', 'b'):
        print(f"  {bin_labels[bin_key]}...")
        geodesic_polygons[bin_key] = build_geodesic_polygons(
            geodesic_file_dict[bin_key], land_polygon)
        print(f"  {len(geodesic_polygons[bin_key]):,} ocean polygons")

if RUN_EASE2:
    print("Building EASE-Grid 2.0 coordinates and bin mappings...")
    ease2_lon, ease2_lat, ease2_x, ease2_y = build_ease2_lonlat()
    ease2_bin_indices = {}
    n_bins_dict = {}
    for bin_key in ('a', 'b'):
        bin_lons, _ = load_geodesic_centres(geodesic_file_dict[bin_key])
        n_bins_dict[bin_key] = len(bin_lons)
        ease2_bin_indices[bin_key] = build_ease2_bin_mapping(
            ease2_lon, ease2_lat, geodesic_file_dict[bin_key])
        print(f"  {bin_labels[bin_key]}: mapped {ease2_bin_indices[bin_key].shape[0]:,} cells")

# Process each data class
for class_label, input_dir in DATA_CLASSES.items():
    print(f"\n{'='*60}")
    print(f"Data class: {class_label}  ({input_dir})")
    print(f"{'='*60}")

    files = sorted(input_dir.glob('*.nc'))
    if not files:
        print("  no .nc files found, skipping")
        continue

    print(f"  collecting records from {len(files)} files...")
    records = collect_records(files, bin_arrays, xyz_grid, grid_indices)
    if not records:
        print("  no valid records, skipping")
        continue
    print(f"  {len(records):,} total profiles")

    class_plots_root = PLOTS_ROOT / class_label
    start_year = detect_start_year(records, 'a')
    end_year = 2025
    print(f"  start year: {start_year}, end year: {end_year}")

    for bin_key, bin_label in bin_labels.items():
        title_base = f'{class_label} — probability of >= N profiles per bin per month'

        # --- All-years combined ---
        probs_all = compute_visit_probabilities(records, bin_key)

        if RUN_POLYGON:
            plot_polygon_maps(
                probs_all, geodesic_polygons[bin_key], bin_label, class_label,
                title=f'{title_base} (all years)',
                out_dir=class_plots_root / 'bin_visit_probability_maps',
                file_stem=f'prob_map_{bin_key}_all_years',
            )

        if RUN_EASE2:
            grids_all = probs_to_ease2_grid(probs_all, ease2_bin_indices[bin_key], n_bins_dict[bin_key])
            plot_ease2_maps(
                grids_all, ease2_x, ease2_y, bin_label, class_label,
                title=f'{title_base} (all years)',
                out_dir=class_plots_root / 'bin_visit_probability_maps_ease2',
                file_stem=f'prob_map_ease2_{bin_key}_all_years',
            )

        # --- Per-year ---
        if start_year is None:
            continue
        for year in range(start_year, end_year + 1):
            probs_yr = compute_visit_probabilities_year(records, bin_key, year)
            if not probs_yr:
                continue

            if RUN_POLYGON:
                plot_polygon_maps(
                    probs_yr, geodesic_polygons[bin_key], bin_label, class_label,
                    title=f'{title_base} ({year})',
                    out_dir=class_plots_root / 'bin_visit_probability_maps' / 'per_year',
                    file_stem=f'prob_map_{bin_key}_{year}',
                )

            if RUN_EASE2:
                grids_yr = probs_to_ease2_grid(probs_yr, ease2_bin_indices[bin_key], n_bins_dict[bin_key])
                plot_ease2_maps(
                    grids_yr, ease2_x, ease2_y, bin_label, class_label,
                    title=f'{title_base} ({year})',
                    out_dir=class_plots_root / 'bin_visit_probability_maps_ease2' / 'per_year',
                    file_stem=f'prob_map_ease2_{bin_key}_{year}',
                )

print("\nDone.")
