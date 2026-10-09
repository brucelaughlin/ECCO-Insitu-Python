"""
Diagnostic plots for spatial density scaling (step 11).

Reproduces the figures from the original MATLAB script
update_weights_based_on_spatial_density_TRYTOIMPLEMENT.m, but per observation
type: step 11 counts and scales T and S separately, so every figure is made
once per variable in SCALED_VARS (filenames/titles tagged _T, _S). A profile
counts for V if it has any nonzero prof_Vweight — same rule as step 11.

  1. Profile count in 1°×1° lat/lon boxes  (raw + log10)
  2. Profile count in 10242 geodesic bins   (global, NH, SH)
  3. Probability of ≥1 profile per month   (global, NH, SH)
  4. Mean scaling factor                   (global, NH, SH)
     — two variants: including and excluding months with no profiles

All figures are saved as PNG under <processed_root>/figures/figures_<ts>/.

Run standalone:
    python plot_spatial_scaling.py /path/to/run_dir
    python plot_spatial_scaling.py /path/to/run_dir --window 10day
    python plot_spatial_scaling.py /path/to/run_dir --dpi 150
"""

import argparse
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature

# Reuse count/window logic from step11 — no duplication
from step11_spatial_scaling import (
    _count_profiles, _has_obs, _WINDOW_FN, N_BINS, SCALED_VARS
)

# ---------------------------------------------------------------------------
# Geodesic bin centres
# ---------------------------------------------------------------------------

_GEO_CSV = Path('/Users/brucel/ecco/yip/sample_data/ecco-insitu/sweet_gdrive/'
                'geodesic/10242_bin_locations.csv')


def _load_bin_centres():
    data = np.genfromtxt(_GEO_CSV, delimiter=',')
    lons = data[:, 0]   # column 0 = lon
    lats = data[:, 1]   # column 1 = lat
    # bin IDs are 1-based; array index = bin_id - 1
    return lons, lats


# ---------------------------------------------------------------------------
# Derive count_per_bin and scaling arrays from raw counts dict
# ---------------------------------------------------------------------------

def _build_arrays(counts, nc_files):
    """
    From step 11's per-variable counts {var: {(bin_id, window_key): n}},
    plus a re-read of lat/lon for the 1°×1° map, build per variable:
      count_per_bin  : (N_BINS,) total profiles with V obs, ever
      count_by_window: dict {bin_id: {window_key: n}}
      lats, lons     : positions of profiles with V obs (for 1deg map)
    Returns {var: (lats, lons, count_per_bin, count_by_window)}.
    """
    import xarray as xr
    pos = {v: ([], []) for v in SCALED_VARS}
    for f in nc_files:
        try:
            ds = xr.open_dataset(f, mask_and_scale=False)
        except Exception:
            continue
        if not {'prof_lat', 'prof_lon'}.issubset(ds.data_vars):
            ds.close()
            continue
        lats = np.asarray(ds['prof_lat'], dtype=np.float32)
        lons = np.asarray(ds['prof_lon'], dtype=np.float32)
        for v in SCALED_VARS:
            if f'prof_{v}weight' in ds:
                has = _has_obs(ds[f'prof_{v}weight'])
                pos[v][0].append(lats[has])
                pos[v][1].append(lons[has])
        ds.close()

    out = {}
    for v in SCALED_VARS:
        count_by_window = {}   # {bin_id (1-based): {window_key: count}}
        for (bid, wk), n in counts[v].items():
            count_by_window.setdefault(bid, {})[wk] = n

        count_per_bin = np.zeros(N_BINS, dtype=np.float64)
        for bid, wdict in count_by_window.items():
            count_per_bin[bid - 1] = sum(wdict.values())

        lats = np.concatenate(pos[v][0]) if pos[v][0] else np.array([])
        lons = np.concatenate(pos[v][1]) if pos[v][1] else np.array([])
        out[v] = (lats, lons, count_per_bin, count_by_window)
    return out


def _build_scaling_arrays(count_by_window, n_windows_total):
    """
    Returns:
      sf_mean_with_empty   : (N_BINS,) mean scaling factor, empty windows = 1.0
      sf_mean_no_empty     : (N_BINS,) mean scaling factor, empty windows excluded
      prob_1plus           : (N_BINS,) fraction of windows with ≥1 profile
    """
    sf_mean_with = np.zeros(N_BINS, dtype=np.float64)
    sf_mean_no   = np.full(N_BINS, np.nan, dtype=np.float64)
    prob         = np.zeros(N_BINS, dtype=np.float64)

    for bid, wdict in count_by_window.items():
        idx = bid - 1
        n_occupied = len(wdict)
        sf_occupied = sum(1.0 / c for c in wdict.values())

        # with empty: sum of (1/c for occupied) + 1.0 * n_empty, divided by total
        n_empty = n_windows_total - n_occupied
        sf_mean_with[idx] = (sf_occupied + float(n_empty)) / n_windows_total

        # without empty: mean only over occupied windows
        sf_mean_no[idx] = sf_occupied / n_occupied if n_occupied > 0 else np.nan

        prob[idx] = n_occupied / n_windows_total

    return sf_mean_with, sf_mean_no, prob


# ---------------------------------------------------------------------------
# Plot helpers
# ---------------------------------------------------------------------------

_CMAP = plt.cm.jet


_MARKER_SIZE = 1.2

def _scatter_geo(ax, lons, lats, values, cmap, vmin, vmax, title, s=_MARKER_SIZE):
    ax.set_global()
    ax.add_feature(cfeature.LAND, facecolor='#dddddd', zorder=0)
    ax.coastlines(linewidth=0.5, zorder=1)
    values = np.asarray(values)
    mask = np.isfinite(values) & (values > 0)
    sc = ax.scatter(lons[mask], lats[mask], c=values[mask], s=s, cmap=cmap,
                    vmin=vmin, vmax=vmax,
                    transform=ccrs.PlateCarree(), zorder=2, rasterized=True)
    ax.set_title(title, fontsize=10)
    return sc


def _make_fig(projection, title, figsize=(12, 6)):
    fig = plt.figure(figsize=figsize)
    ax  = fig.add_subplot(1, 1, 1, projection=projection)
    fig.suptitle(title, fontsize=11)
    return fig, ax


def _save(fig, path, dpi):
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=dpi, bbox_inches='tight')
    plt.close(fig)
    print(f"  saved: {path.name}")


# ---------------------------------------------------------------------------
# Individual plot functions
# ---------------------------------------------------------------------------

def plot_1deg_count(var, all_lats, all_lons, fig_dir, run_label, dpi):
    """Figure 1: profile count in 1°×1° lat/lon boxes."""
    lat_edges = np.arange(-90, 91, 1)
    lon_edges = np.arange(-180, 181, 1)
    counts, _, _ = np.histogram2d(all_lats, all_lons,
                                   bins=[lat_edges, lon_edges])

    lons_c = 0.5 * (lon_edges[:-1] + lon_edges[1:])
    lats_c = 0.5 * (lat_edges[:-1] + lat_edges[1:])
    lon_g, lat_g = np.meshgrid(lons_c, lats_c)

    mask   = counts > 0
    vmax   = 2 * np.median(counts[mask]) if mask.any() else 1

    fig, axes = plt.subplots(1, 2, figsize=(16, 5),
                              subplot_kw={'projection': ccrs.PlateCarree()})
    for ax in axes:
        ax.set_global()
        ax.add_feature(cfeature.LAND, facecolor='#dddddd', zorder=0)
        ax.coastlines(linewidth=0.4, zorder=1)

    pcm0 = axes[0].pcolormesh(lon_g, lat_g, counts, cmap=_CMAP,
                               vmin=0, vmax=vmax,
                               transform=ccrs.PlateCarree(), rasterized=True)
    axes[0].set_title(f'{var} profile count (1°×1° boxes)', fontsize=10)
    plt.colorbar(pcm0, ax=axes[0], orientation='horizontal', pad=0.04, shrink=0.8)

    log_counts = np.where(mask, np.log10(counts), 0)
    pcm1 = axes[1].pcolormesh(lon_g, lat_g, log_counts, cmap=_CMAP,
                               transform=ccrs.PlateCarree(), rasterized=True)
    axes[1].set_title(f'log₁₀ {var} profile count (1°×1° boxes)', fontsize=10)
    plt.colorbar(pcm1, ax=axes[1], orientation='horizontal', pad=0.04, shrink=0.8)

    fig.suptitle(f'{var} profile count — {run_label}', fontsize=11)
    _save(fig, fig_dir / f'profile_count_1deg_{var}_{run_label}.png', dpi)


def _bin_scatter_panel(lons_g, lats_g, values, cmap, vmin, vmax,
                        title, projection, fig_dir, fname, run_label, dpi):
    fig, ax = _make_fig(projection, f'{title} — {run_label}')
    sc = _scatter_geo(ax, lons_g, lats_g, values, cmap, vmin, vmax, title)
    plt.colorbar(sc, ax=ax, orientation='horizontal', pad=0.04, shrink=0.8)
    _save(fig, fig_dir / fname, dpi)


def plot_bin_counts(var, lons_g, lats_g, count_per_bin, fig_dir, run_label, dpi):
    """Figures 2a-c: profile count in geodesic bins — global, NH, SH."""
    valid  = count_per_bin[count_per_bin > 0]
    vmax   = 10 * np.median(valid) if len(valid) else 1

    for proj, suffix, title_sfx in [
        (ccrs.Robinson(),       'global', 'Global'),
        (ccrs.Orthographic(0,  90), 'NH',     'Northern Hemisphere'),
        (ccrs.Orthographic(0, -90), 'SH',     'Southern Hemisphere'),
    ]:
        _bin_scatter_panel(
            lons_g, lats_g, count_per_bin, _CMAP, 0, vmax,
            f'{var} profile count in 10242 bins — {title_sfx}',
            proj, fig_dir,
            f'profile_count_10242bins_{var}_{suffix}_{run_label}.png',
            run_label, dpi,
        )

    # log10 global
    log_count = np.where(count_per_bin > 0, np.log10(count_per_bin), 0)
    _bin_scatter_panel(
        lons_g, lats_g, log_count, _CMAP, 1,
        min(4, int(np.ceil(log_count.max())) if log_count.max() > 0 else 4),
        f'log₁₀ {var} profile count in 10242 bins',
        ccrs.Robinson(), fig_dir,
        f'profile_count_10242bins_log10_{var}_{run_label}.png',
        run_label, dpi,
    )


def plot_probability(var, lons_g, lats_g, prob, fig_dir, run_label, dpi):
    """Figures 3a-c: P(≥1 profile per window)."""
    for proj, suffix, title_sfx in [
        (ccrs.Robinson(),       'global', 'Global'),
        (ccrs.Orthographic(0,  90), 'NH',     'Northern Hemisphere'),
        (ccrs.Orthographic(0, -90), 'SH',     'Southern Hemisphere'),
    ]:
        _bin_scatter_panel(
            lons_g, lats_g, prob, plt.cm.viridis, 0, 1,
            f'P(≥1 {var} profile / window) — {title_sfx}',
            proj, fig_dir,
            f'prob_1plus_{var}_{suffix}_{run_label}.png',
            run_label, dpi,
        )


def plot_mean_scaling(var, lons_g, lats_g, sf_with, sf_no, fig_dir, run_label, dpi):
    """Figures 4a-b: mean scaling factor (with / without empty windows)."""
    for values, label, tag in [
        (sf_with, 'mean scaling factor (empty windows = 1)', 'with_empty'),
        (np.where(np.isnan(sf_no), 0, sf_no), 'mean scaling factor (occupied windows only)', 'no_empty'),
    ]:
        for proj, suffix, title_sfx in [
            (ccrs.Robinson(),       'global', 'Global'),
            (ccrs.Orthographic(0,  90), 'NH',     'NH'),
            (ccrs.Orthographic(0, -90), 'SH',     'SH'),
        ]:
            _bin_scatter_panel(
                lons_g, lats_g, values, _CMAP, 0, 1,
                f'{var} {label} — {title_sfx}',
                proj, fig_dir,
                f'mean_scaling_factor_{var}_{tag}_{suffix}_{run_label}.png',
                run_label, dpi,
            )


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main(run_dir, window_mode='monthly', dpi=150):
    run_dir = Path(run_dir)
    if not run_dir.is_dir():
        raise ValueError(f"run_dir does not exist: {run_dir}")

    window_fn = _WINDOW_FN.get(window_mode)
    if window_fn is None:
        raise ValueError(f"Unknown window_mode '{window_mode}'. Choose from: {list(_WINDOW_FN)}")

    nc_files = sorted(run_dir.rglob('*.nc'))
    if not nc_files:
        print(f"No .nc files found under {run_dir}")
        return

    run_label  = run_dir.name

    # Extract timestamp from run directory name (e.g. ..._20261007_144529)
    import re as _re
    _ts_match = _re.search(r'_(\d{8}_\d{6})$', run_dir.name)
    _ts = _ts_match.group(1) if _ts_match else run_dir.name
    fig_dir = run_dir.parent / 'figures' / f'figures_{_ts}'

    print(f"\n{'='*60}")
    print(f"Spatial scaling diagnostic plots")
    print(f"  run_dir : {run_dir}")
    print(f"  window  : {window_mode}")
    print(f"  files   : {len(nc_files)}")
    print(f"  fig_dir : {fig_dir}")
    print(f"{'='*60}\n")

    print("Loading bin centres...")
    lons_g, lats_g = _load_bin_centres()

    print("Counting profiles per variable (step 11 pass 1)...")
    counts = _count_profiles(nc_files, window_fn)

    print("Building arrays...")
    per_var = _build_arrays(counts, nc_files)

    # Number of distinct windows in the dataset (any variable) — shared
    # denominator so T and S probabilities are comparable
    all_windows = set()
    for v in SCALED_VARS:
        all_windows.update(wk for (_, wk) in counts[v])
    n_windows_total = len(all_windows)
    print(f"  {n_windows_total} distinct ({window_mode}) windows in dataset")

    for v in SCALED_VARS:
        all_lats, all_lons, count_per_bin, count_by_window = per_var[v]
        print(f"\nGenerating {v} figures ({int(count_per_bin.sum()):,} profiles)...")
        if not count_by_window:
            print(f"  no {v} observations, skipping")
            continue

        sf_with, sf_no, prob = _build_scaling_arrays(count_by_window, n_windows_total)

        if len(all_lats):
            print("  1deg count map...")
            plot_1deg_count(v, all_lats, all_lons, fig_dir, run_label, dpi)

        print("  Geodesic bin count maps...")
        plot_bin_counts(v, lons_g, lats_g, count_per_bin, fig_dir, run_label, dpi)

        print("  Probability maps...")
        plot_probability(v, lons_g, lats_g, prob, fig_dir, run_label, dpi)

        print("  Mean scaling factor maps...")
        plot_mean_scaling(v, lons_g, lats_g, sf_with, sf_no, fig_dir, run_label, dpi)

    print(f"\nDone. {len(list(fig_dir.glob('*.png')))} figures saved to {fig_dir}\n")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Diagnostic plots for step 11 spatial density scaling.')
    parser.add_argument('run_dir',
                        help='Root directory of a completed NCEI run')
    parser.add_argument('--window', default='monthly', choices=list(_WINDOW_FN),
                        help='Window mode (default: monthly)')
    parser.add_argument('--dpi', type=int, default=150,
                        help='Figure resolution in DPI (default: 150)')
    args = parser.parse_args()
    main(args.run_dir, window_mode=args.window, dpi=args.dpi)
