"""
For each month, count the number of bins visited by at least 1, at least 2,
or at least 3 profiles. Plot as line plots over time (one point per month).

Separate plot sets are produced for each bin resolution (a/b) and each
profile type (all combined, A, D, R).
"""

import sys
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from pathlib import Path
from collections import defaultdict
from scipy.interpolate import griddata

sys.path.append('/Users/brucel/ecco/yip/ECCO-Insitu-Python')
import tools

# ==============================================================================
# Switches
# ==============================================================================

RUN_RAW       = True
RUN_PROCESSED = False  # set True once NCEI pipeline output is available

# ==============================================================================
# Paths
# ==============================================================================

raw_pfl_dir   = Path('/Users/brucel/ecco/yip/profile_data/Interp_Profiles/PFL')
processed_dir = Path('/Users/brucel/ecco/yip/profile_files_NCEI_processed/PFL')

sphere_bin_dir = '/Users/brucel/ecco/yip/sample_data/ecco-insitu/sweet_gdrive/grid_llc90/sphere_point_distribution'
grid_dir       = '/Users/brucel/ecco/yip/sample_data/ecco-insitu/sweet_gdrive/grid_llc90'

output_dir = Path('/Users/brucel/ecco/yip/profile_data/z_profile_file_analysis/plots/bin_visit_threshold_timeseries')
output_dir.mkdir(parents=True, exist_ok=True)

bin_labels = {'a': '10242 bins', 'b': '02562 bins'}

# ==============================================================================
# Grid / bin helpers (shared with plot_bin_occupancy_per_month.py)
# ==============================================================================

def load_bin_arrays(sphere_bin_dir):
    import os
    import numpy as np
    tile_shape = [90, 13*90, 1]
    mform = '>f4'
    bin_arrays = {}
    for key, fname in [('a', 'llc090_sphere_point_n_10242_ids.bin'),
                       ('b', 'llc090_sphere_point_n_02562_ids.bin')]:
        fpath = os.path.join(sphere_bin_dir, fname)
        with open(fpath, 'rb') as fid:
            arr = np.fromfile(fid, dtype=mform)
        arr = arr.reshape((tile_shape[0], tile_shape[1], tile_shape[2]), order='F')
        bin_arrays[key] = arr.ravel()
    return bin_arrays


def build_llc90_griddata_inputs(grid_dir):
    lon_90, lat_90, bathy_90, X_90, Y_90, Z_90 = tools.load_llc90_grid_step2(grid_dir)
    xyz_grid = np.column_stack((X_90.ravel(), Y_90.ravel(), Z_90.ravel()))
    indices = np.arange(X_90.size)
    return xyz_grid, indices


def assign_bin_ids_from_latlon(prof_lon, prof_lat, xyz_grid, grid_indices, bin_arrays):
    deg2rad = np.pi / 180.0
    valid_mask = np.asarray(
        tools.sph2cart_returnValidMaskOnly(prof_lon * deg2rad, prof_lat * deg2rad, 1),
        dtype=bool)
    x, y, z = tools.sph2cart(prof_lon[valid_mask] * deg2rad,
                              prof_lat[valid_mask] * deg2rad, 1)
    xyz_prof = np.column_stack((np.asarray(x), np.asarray(y), np.asarray(z)))
    llc_indices = griddata(xyz_grid, grid_indices, xyz_prof, method='nearest').astype(int)
    result = {}
    for key, bin_arr in bin_arrays.items():
        ids = np.full(len(prof_lon), np.nan)
        ids[valid_mask] = bin_arr[llc_indices]
        result[key] = ids
    return result

# ==============================================================================
# Data collection
# ==============================================================================

def collect_records(file_list, bin_arrays, xyz_grid=None, grid_indices=None,
                    use_precomputed_bins=False):
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

        prof_type = fpath.stem.split('_')[-1][0]
        yyyymmdd = ds['prof_YYYYMMDD'].values.astype(int)
        months = (yyyymmdd % 10000) // 100
        years  = yyyymmdd // 10000

        if use_precomputed_bins:
            bin_id_a = ds['prof_bin_id_a'].values
            bin_id_b = ds['prof_bin_id_b'].values
        else:
            prof_lon = ds['prof_lon'].values.astype(float)
            prof_lat = ds['prof_lat'].values.astype(float)
            assigned = assign_bin_ids_from_latlon(
                prof_lon, prof_lat, xyz_grid, grid_indices, bin_arrays)
            bin_id_a = assigned['a']
            bin_id_b = assigned['b']

        for i in range(len(months)):
            records.append({
                'year':      int(years[i]),
                'month':     int(months[i]),
                'bin_id_a':  float(bin_id_a[i]),
                'bin_id_b':  float(bin_id_b[i]),
                'prof_type': prof_type,
            })
        ds.close()
        print(f"  loaded {fpath.name}: {len(months)} profiles")
    return records

# ==============================================================================
# Counting
# ==============================================================================

def count_threshold_timeseries(records, bin_key):
    """
    For each (year, month), count bins visited by >= 1, >= 2, >= 3 profiles.
    Returns dict: (year, month) -> {'1+': N, '2+': N, '3+': N}
    """
    ym_bin_counts = defaultdict(lambda: defaultdict(int))
    for r in records:
        bid = r[f'bin_id_{bin_key}']
        if np.isnan(bid):
            continue
        ym_bin_counts[(r['year'], r['month'])][bid] += 1

    result = {}
    for ym, bin_count_dict in ym_bin_counts.items():
        counts = list(bin_count_dict.values())
        result[ym] = {
            '1+': sum(1 for c in counts if c >= 1),
            '2+': sum(1 for c in counts if c >= 2),
            '3+': sum(1 for c in counts if c >= 3),
        }
    return result

# ==============================================================================
# Plotting
# ==============================================================================

def plot_threshold_timeseries(ts_dict, title, output_path):
    """
    Line plot: one point per month, x = fractional year, y = bin count.
    Three lines per bin resolution: >= 1, >= 2, >= 3 profiles per bin.
    Two panels (one per bin resolution), stacked vertically.
    """
    categories = ['1+', '2+', '3+']
    colors     = ['#4878d0', '#ee854a', '#6acc65']
    fmt = ticker.FuncFormatter(lambda v, _: f'{v/1000:.1f}k' if v >= 1000 else str(int(v)))

    fig, axes = plt.subplots(2, 1, figsize=(14, 7), sharex=True)

    for ax, (bin_key, bin_label) in zip(axes, bin_labels.items()):
        ts = ts_dict[bin_key]
        sorted_keys = sorted(ts.keys())
        x = [ym[0] + (ym[1] - 1) / 12.0 for ym in sorted_keys]

        for cat, color in zip(categories, colors):
            y = [ts[ym][cat] for ym in sorted_keys]
            ax.plot(x, y, linewidth=1.2, color=color, label=f'>= {cat[0]} profile(s)/bin/month')

        ax.set_ylabel('Number of bins')
        ax.set_title(bin_label, fontsize=9)
        ax.legend(fontsize=8)
        ax.yaxis.set_major_formatter(fmt)
        ax.grid(axis='y', linewidth=0.4, alpha=0.5)

    axes[-1].set_xlabel('Year')
    fig.suptitle(title, fontsize=11)
    fig.tight_layout()
    fig.savefig(output_path, dpi=150)
    plt.close(fig)
    print(f"saved: {output_path}")

# ==============================================================================
# Main
# ==============================================================================

print("Loading bin arrays...")
bin_arrays = load_bin_arrays(sphere_bin_dir)

datasets = []

if RUN_RAW:
    # PFL_BGC excluded — potential systematic issues would distort the raw count
    print("\nCollecting raw profiles (PFL only)...")
    xyz_grid, grid_indices = build_llc90_griddata_inputs(grid_dir)
    raw_files = sorted(raw_pfl_dir.rglob('*.nc'))
    raw_records = collect_records(raw_files, bin_arrays, xyz_grid=xyz_grid,
                                  grid_indices=grid_indices, use_precomputed_bins=False)
    print(f"  total raw profiles: {len(raw_records):,}")
    datasets.append(('raw', raw_records))

if RUN_PROCESSED:
    print("\nCollecting processed profiles (PFL)...")
    proc_files = sorted(processed_dir.rglob('*.nc'))
    proc_records = collect_records(proc_files, bin_arrays, use_precomputed_bins=True)
    print(f"  total processed profiles: {len(proc_records):,}")
    datasets.append(('processed', proc_records))

for label, records in datasets:
    record_sets = [('all', records)] + [
        (f'type{t}', [r for r in records if r.get('prof_type') == t])
        for t in ('A', 'D', 'R')
        if any(r.get('prof_type') == t for r in records)
    ]

    for type_label, recs in record_sets:
        print(f"\nCounting threshold occupancy ({label}, {type_label}): {len(recs):,} profiles")
        ts_dict = {}
        for bin_key in ('a', 'b'):
            ts_dict[bin_key] = count_threshold_timeseries(recs, bin_key)

        plot_threshold_timeseries(
            ts_dict,
            title=f'Argo PFL — bins with >= N profiles per month ({label}, {type_label})',
            output_path=output_dir / f'argo_bin_threshold_timeseries_{label}_{type_label}.png',
        )

print("\nDone.")
