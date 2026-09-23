"""
Plot the number of geodesic bins occupied by at least one Argo profile per month.

Two bin resolutions: 10242 (prof_bin_id_a) and 02562 (prof_bin_id_b).
Two data sources: raw (before NCEI chain) and processed (after NCEI chain).
Two time views: calendar month (Jan-Dec, aggregated across all years) and
                month-by-year time series.

Output: 4 figures saved to ../plots/data_coverage/
"""

import os
import sys
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from pathlib import Path
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

sphere_bin_dir  = '/Users/brucel/ecco/yip/sample_data/ecco-insitu/sweet_gdrive/grid_llc90/sphere_point_distribution'
grid_dir        = '/Users/brucel/ecco/yip/sample_data/ecco-insitu/sweet_gdrive/grid_llc90'

output_dir = Path('/Users/brucel/ecco/yip/profile_data/z_profile_file_analysis/plots/bin_occupancy')
output_dir.mkdir(parents=True, exist_ok=True)

bin_labels = {'a': '10242 bins', 'b': '02562 bins'}

# ==============================================================================
# Load geodesic bin ID arrays (LLC90 grid point → bin ID)
# ==============================================================================

def load_bin_arrays(sphere_bin_dir):
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
    """Return xyz_grid and flattened_monotonic_grid_indices for LLC90."""
    deg2rad = np.pi / 180.0
    lon_90, lat_90, bathy_90, X_90, Y_90, Z_90 = tools.load_llc90_grid_step2(grid_dir)
    xyz_grid = np.column_stack((X_90.ravel(), Y_90.ravel(), Z_90.ravel()))
    indices = np.arange(X_90.size)
    return xyz_grid, indices


def assign_bin_ids_from_latlon(prof_lon, prof_lat, xyz_grid, grid_indices, bin_arrays):
    """Assign bin IDs to profiles given their lat/lon."""
    deg2rad = np.pi / 180.0
    valid_mask = tools.sph2cart_returnValidMaskOnly(
        prof_lon * deg2rad, prof_lat * deg2rad, 1)
    valid_mask = np.asarray(valid_mask, dtype=bool)

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
# Collect data from files
# ==============================================================================

def collect_records(file_list, bin_arrays, xyz_grid=None, grid_indices=None,
                    use_precomputed_bins=False):
    """
    Returns a list of dicts, one per profile, with keys:
      'year', 'month', 'bin_id_a', 'bin_id_b'
    """
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

        # Extract profile type (A, D, or R) from filename, e.g. ARGO_WO_2005_PFL_D.nc
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


def count_occupied_bins(records, bin_key):
    """
    Returns two dicts:
      calendar_counts[month]       = number of unique bins occupied in that calendar month
      timeseries_counts[(year,month)] = number of unique bins occupied
    """
    from collections import defaultdict
    calendar = defaultdict(set)
    timeseries = defaultdict(set)

    for r in records:
        bid = r[f'bin_id_{bin_key}']
        if np.isnan(bid):
            continue
        m = r['month']
        ym = (r['year'], r['month'])
        calendar[m].add(bid)
        timeseries[ym].add(bid)

    calendar_counts   = {m: len(s) for m, s in calendar.items()}
    timeseries_counts = {ym: len(s) for ym, s in timeseries.items()}
    return calendar_counts, timeseries_counts


def count_visit_frequency(records, bin_key):
    """
    For each (year, month), count how many bins were visited exactly 1x, 2x, 3+ times.

    Returns two dicts keyed by month and (year,month):
      calendar_freq[month]    = {'1x': N, '2x': N, '3+': N}  (summed across all years)
      timeseries_freq[(y,m)]  = {'1x': N, '2x': N, '3+': N}
    """
    from collections import defaultdict
    # bin_id -> count, per (year, month)
    ym_bin_counts = defaultdict(lambda: defaultdict(int))

    for r in records:
        bid = r[f'bin_id_{bin_key}']
        if np.isnan(bid):
            continue
        ym_bin_counts[(r['year'], r['month'])][bid] += 1

    timeseries_freq = {}
    for ym, bin_count_dict in ym_bin_counts.items():
        counts = list(bin_count_dict.values())
        timeseries_freq[ym] = {
            '1x':  sum(1 for c in counts if c == 1),
            '2x':  sum(1 for c in counts if c == 2),
            '3+':  sum(1 for c in counts if c >= 3),
        }

    # Aggregate across years for calendar view
    calendar_freq = {m: {'1x': 0, '2x': 0, '3+': 0} for m in range(1, 13)}
    for (y, m), freq in timeseries_freq.items():
        for cat in ('1x', '2x', '3+'):
            calendar_freq[m][cat] += freq[cat]

    return calendar_freq, timeseries_freq

# ==============================================================================
# Plotting
# ==============================================================================

MONTH_NAMES = ['Jan','Feb','Mar','Apr','May','Jun',
               'Jul','Aug','Sep','Oct','Nov','Dec']


def plot_calendar(calendar_counts_dict, title, output_path):
    """Bar chart: x = calendar month, y = occupied bins. One line per bin resolution."""
    fig, ax = plt.subplots(figsize=(9, 4))
    x = np.arange(1, 13)
    width = 0.35
    offsets = [-width/2, width/2]

    for offset, (bin_key, label) in zip(offsets, bin_labels.items()):
        counts = calendar_counts_dict[bin_key]
        y = [counts.get(m, 0) for m in range(1, 13)]
        ax.bar(x + offset, y, width, label=label)

    ax.set_xticks(x)
    ax.set_xticklabels(MONTH_NAMES)
    ax.set_xlabel('Month')
    ax.set_ylabel('Occupied bins')
    ax.set_title(title)
    ax.legend()
    ax.yaxis.set_major_formatter(ticker.FuncFormatter(lambda v, _: f'{int(v):,}'))
    fig.tight_layout()
    fig.savefig(output_path, dpi=150)
    plt.close(fig)
    print(f"saved: {output_path}")


def plot_timeseries(timeseries_counts_dict, title, output_path):
    """Line plot: x = year-month, y = occupied bins. One line per bin resolution."""
    fig, ax = plt.subplots(figsize=(14, 4))

    for bin_key, label in bin_labels.items():
        counts = timeseries_counts_dict[bin_key]
        sorted_keys = sorted(counts.keys())
        # fractional year for x axis
        x = [ym[0] + (ym[1] - 1) / 12.0 for ym in sorted_keys]
        y = [counts[ym] for ym in sorted_keys]
        ax.plot(x, y, linewidth=0.8, label=label)

    ax.set_xlabel('Year')
    ax.set_ylabel('Occupied bins')
    ax.set_title(title)
    ax.legend()
    ax.yaxis.set_major_formatter(ticker.FuncFormatter(lambda v, _: f'{int(v):,}'))
    fig.tight_layout()
    fig.savefig(output_path, dpi=150)
    plt.close(fig)
    print(f"saved: {output_path}")


def _render_year_grid_page(years_page, timeseries_freq, bin_key, bin_label,
                           year_to_regime, year_to_ymax, title, output_path,
                           page_num=None, n_pages=None):
    """Render one page of the year grid: grouped 1x/2x/3+ bars per subplot."""
    ncols = 4
    nrows = int(np.ceil(len(years_page) / ncols))
    categories = ['1x', '2x', '3+']
    colors = ['#4878d0', '#ee854a', '#6acc65']
    width = 0.25
    x = np.arange(1, 13)
    offsets = np.array([-1, 0, 1]) * width
    fmt = ticker.FuncFormatter(lambda v, _: f'{v/1000:.1f}k' if v >= 1000 else str(int(v)))

    fig, axes = plt.subplots(nrows, ncols, figsize=(ncols * 3.5, nrows * 2.4),
                             sharey=False, sharex=False)
    axes_flat = np.array(axes).ravel()

    for ax_idx, year in enumerate(years_page):
        ax = axes_flat[ax_idx]
        for cat, color, offset in zip(categories, colors, offsets):
            y = np.array([timeseries_freq.get((year, m), {}).get(cat, 0)
                          for m in range(1, 13)], dtype=float)
            ax.bar(x + offset, y, width, color=color, label=cat)

        ymax = year_to_ymax[year]
        regime = year_to_regime[year]
        ax.set_ylim(0, ymax * 1.1)
        ax.set_title(str(year), fontsize=8)
        ax.set_xticks(x)
        ax.set_xticklabels(MONTH_NAMES, fontsize=5, rotation=45, ha='right')
        ax.tick_params(axis='y', labelsize=6)
        ax.yaxis.set_major_formatter(fmt)
        ax.text(0.02, 0.97, f'regime {regime}: y-max {fmt(ymax, None)}',
                transform=ax.transAxes, fontsize=6, va='top', color='dimgray')

    for ax_idx in range(len(years_page), len(axes_flat)):
        axes_flat[ax_idx].set_visible(False)

    handles, labels_legend = axes_flat[0].get_legend_handles_labels()
    fig.legend(handles, labels_legend, title='profiles/bin/month', loc='lower right', fontsize=8)
    page_str = f' — page {page_num}/{n_pages}' if page_num is not None and n_pages is not None else ''
    fig.suptitle(f'{title}\n{bin_label}{page_str}', fontsize=10)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    fig.savefig(output_path, dpi=150)
    plt.close(fig)
    print(f"saved: {output_path}")


def plot_visit_freq_by_year_grid(timeseries_freq_dict, label, title_base, output_path_base):
    """
    For each bin resolution, produce 2 pages of 4x4 subplots (one per year).
    Years are split into up to 3 regimes by tercile of their annual max occupancy.
    Each regime shares a y-scale equal to its group's maximum, so every regime has
    at least one year with bars near the top. Each subplot labels its regime and y-max.
    """
    for bin_key, bin_label in bin_labels.items():
        ts_freq = timeseries_freq_dict[bin_key]
        all_years = sorted({ym[0] for ym in ts_freq})

        def year_max(year):
            return max((sum(ts_freq.get((year, m), {}).values()) for m in range(1, 13)),
                       default=0)

        # Peak bar height = tallest individual category bar across all months.
        def year_peak_bar(year):
            return max(
                (ts_freq.get((year, m), {}).get(cat, 0)
                 for m in range(1, 13) for cat in ('1x', '2x', '3+')),
                default=0
            )

        year_maxes = {y: year_peak_bar(y) for y in all_years}
        sorted_by_max = sorted(all_years, key=lambda y: year_maxes[y])

        # Fixed thresholds for regimes 1 and 2; regime 3 uses data max + 100
        regime3_ymax = max(year_maxes.values()) + 100
        regime_bounds = [50, 200]  # upper bound (inclusive) for regimes 1 and 2

        regime_groups = [[], [], []]
        for y in sorted_by_max:
            if year_maxes[y] <= regime_bounds[0]:
                regime_groups[0].append(y)
            elif year_maxes[y] <= regime_bounds[1]:
                regime_groups[1].append(y)
            else:
                regime_groups[2].append(y)
        regime_groups = [g for g in regime_groups if g]

        year_to_regime = {}
        year_to_ymax = {}
        for y in sorted_by_max:
            if year_maxes[y] <= regime_bounds[0]:
                year_to_regime[y] = 1
                year_to_ymax[y] = regime_bounds[0]
            elif year_maxes[y] <= regime_bounds[1]:
                year_to_regime[y] = 2
                year_to_ymax[y] = regime_bounds[1]
            else:
                year_to_regime[y] = 3
                year_to_ymax[y] = regime3_ymax

        ordered_years = [y for group in regime_groups for y in group]
        mid = len(ordered_years) // 2
        pages = [p for p in [ordered_years[:mid], ordered_years[mid:]] if p]
        n_pages = len(pages)

        for page_num, years_page in enumerate(pages, start=1):
            out = Path(str(output_path_base).format(bin_key=bin_key, page=page_num))
            _render_year_grid_page(
                years_page, ts_freq, bin_key, bin_label,
                year_to_regime, year_to_ymax,
                title=f'{title_base} ({label})', output_path=out,
                page_num=page_num, n_pages=n_pages,
            )


def plot_visit_frequency_calendar(calendar_freq, bin_label, title, output_path):
    """Grouped bar chart: x = calendar month, grouped by 1x/2x/3+ visit frequency."""
    fig, ax = plt.subplots(figsize=(9, 4))
    x = np.arange(1, 13)
    categories = ['1x', '2x', '3+']
    colors = ['#4878d0', '#ee854a', '#6acc65']
    width = 0.25
    offsets = np.array([-1, 0, 1]) * width

    for cat, color, offset in zip(categories, colors, offsets):
        y = np.array([calendar_freq[m][cat] for m in range(1, 13)], dtype=float)
        ax.bar(x + offset, y, width, label=cat, color=color)

    ax.set_xticks(x)
    ax.set_xticklabels(MONTH_NAMES)
    ax.set_xlabel('Month')
    ax.set_ylabel('Number of bins')
    ax.set_title(f'{title}\n{bin_label}')
    ax.legend(title='profiles/bin/month')
    ax.yaxis.set_major_formatter(ticker.FuncFormatter(lambda v, _: f'{int(v):,}'))
    fig.tight_layout()
    fig.savefig(output_path, dpi=150)
    plt.close(fig)
    print(f"saved: {output_path}")


def plot_visit_frequency_timeseries(timeseries_freq, bin_label, title, output_path):
    """Stacked area plot: x = year-month, stacked by 1x/2x/3+ visit frequency."""
    fig, ax = plt.subplots(figsize=(14, 4))
    sorted_keys = sorted(timeseries_freq.keys())
    x = [ym[0] + (ym[1] - 1) / 12.0 for ym in sorted_keys]
    categories = ['1x', '2x', '3+']
    colors = ['#4878d0', '#ee854a', '#6acc65']
    ys = [np.array([timeseries_freq[ym][cat] for ym in sorted_keys], dtype=float)
          for cat in categories]

    ax.stackplot(x, ys, labels=categories, colors=colors, alpha=0.85)
    ax.set_xlabel('Year')
    ax.set_ylabel('Number of bins')
    ax.set_title(f'{title}\n{bin_label}')
    ax.legend(title='profiles/bin/month', loc='upper left')
    ax.yaxis.set_major_formatter(ticker.FuncFormatter(lambda v, _: f'{int(v):,}'))
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
    # PFL_BGC is excluded here — BGC files may have a systematic issue that caused
    # them to fail the NCEI chain, so including them in the raw count would make the
    # before/after comparison misleading.
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

def run_plots_for_records(records, label):
    """Run the full suite of occupancy and visit-frequency plots for one record set."""
    print(f"\nCounting occupied bins ({label})...")
    cal_dict = {}
    ts_dict  = {}
    for bin_key in ('a', 'b'):
        cal, ts = count_occupied_bins(records, bin_key)
        cal_dict[bin_key] = cal
        ts_dict[bin_key]  = ts

    plot_calendar(
        cal_dict,
        title=f'Argo PFL — occupied geodesic bins per calendar month ({label})',
        output_path=output_dir / f'argo_bin_occupancy_calendar_{label}.png'
    )
    plot_timeseries(
        ts_dict,
        title=f'Argo PFL — occupied geodesic bins per month ({label})',
        output_path=output_dir / f'argo_bin_occupancy_timeseries_{label}.png'
    )

    print(f"\nCounting visit frequency ({label})...")
    ts_freq_dict = {}
    for bin_key, bin_label in bin_labels.items():
        cal_freq, ts_freq = count_visit_frequency(records, bin_key)
        ts_freq_dict[bin_key] = ts_freq
        plot_visit_frequency_calendar(
            cal_freq, bin_label,
            title=f'Argo PFL — bin visit frequency per calendar month ({label})',
            output_path=output_dir / f'argo_bin_visit_freq_calendar_{label}_{bin_key}.png'
        )
        plot_visit_frequency_timeseries(
            ts_freq, bin_label,
            title=f'Argo PFL — bin visit frequency per month ({label})',
            output_path=output_dir / f'argo_bin_visit_freq_timeseries_{label}_{bin_key}.png'
        )

    plot_visit_freq_by_year_grid(
        ts_freq_dict, label,
        title_base='Argo PFL — bin visit frequency by year',
        output_path_base=output_dir / f'argo_bin_visit_freq_by_year_{label}_{{bin_key}}_page{{page}}.png'
    )


# --- Count and plot ---
for label, records in datasets:
    # All profiles combined
    run_plots_for_records(records, label)

    # Per profile type (A = ascending, D = descending, R = real-time)
    for prof_type in ('A', 'D', 'R'):
        type_records = [r for r in records if r.get('prof_type') == prof_type]
        if not type_records:
            print(f"  no records found for type {prof_type} in {label}, skipping")
            continue
        print(f"\n  profile type {prof_type}: {len(type_records):,} profiles")
        run_plots_for_records(type_records, f'{label}_type{prof_type}')

print("\nDone.")
