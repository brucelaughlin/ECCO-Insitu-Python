"""
Preprocessing script for the interactive PFL globe app.

Uses post-NCEI step-10 files from profile_files_NCEI_processed/PFL/.
These have precomputed bin IDs (prof_bin_id_a) and per-profile climatology
(prof_Tclim, prof_Sclim) so no LLC90 griddata step is needed.

Per bin, per (year, month), computes nanmean of:
  T, S          — observed profiles
  Tclim, Sclim  — climatology interpolated to each profile's location/time

The app then offers:
  T, S                — bin-averaged observations
  Tclim, Sclim        — bin-averaged climatology
  T−clim, S−clim      — obs minus spatiotemporal climatology (boss request)
  T_anom, S_anom      — obs minus bin's own long-run time mean

Output pickle:
  'depth'        : np.ndarray (97,)
  'bins'         : dict[bin_id -> {
      'bin_id'   : int
      'months'   : dict[(year,month) -> {
          'T','S','Tclim','Sclim': np.ndarray (97,)
          'sources': list of {file, prof_idx, date, lon, lat, type}
      }]
      'T_mean'   : np.ndarray (97,)   long-run time-mean of monthly T
      'S_mean'   : np.ndarray (97,)   long-run time-mean of monthly S
      'prob_1plus': float
  }]
  'bin_centres'  : dict[bin_id -> (lon, lat)]
"""

import pickle
import numpy as np
import pandas as pd
import xarray as xr
from pathlib import Path
from collections import defaultdict

# ==============================================================================
# Paths
# ==============================================================================

PFL_DIR       = Path('/Users/brucel/ecco/yip/profile_files_NCEI_processed/PFL')
GEODESIC_FILE = '/Users/brucel/ecco/yip/sample_data/ecco-insitu/sweet_gdrive/geodesic/10242_bin_locations.csv'
OUTPUT_PICKLE = Path('/Users/brucel/ecco/yip/profile_data/z_profile_file_analysis/pfl_bin_timeseries.pkl')

BIN_ID_VAR  = 'prof_bin_id_a'   # 10242-bin resolution, precomputed by NCEI chain
N_DEPTH_REF = 97                 # full depth grid length

# ==============================================================================
# Helpers
# ==============================================================================

def load_geodesic_centres(geodesic_file):
    df = pd.read_csv(geodesic_file, header=None, names=['lon', 'lat'])
    return {i + 1: (float(df['lon'].iloc[i]), float(df['lat'].iloc[i]))
            for i in range(len(df))}


def pad_to_ref(arr, n_ref=N_DEPTH_REF):
    """Pad a 1-D depth profile to n_ref levels with NaNs at the tail."""
    if len(arr) == n_ref:
        return arr
    out = np.full(n_ref, np.nan)
    out[:len(arr)] = arr
    return out

# ==============================================================================
# Main
# ==============================================================================

print("Loading geodesic bin centres...")
bin_centres = load_geodesic_centres(GEODESIC_FILE)
print(f"  {len(bin_centres):,} bin centres loaded")

# Establish the reference 97-level depth array from the first file that has it
_ref_file = next(
    (f for f in sorted(PFL_DIR.glob('*.nc'))
     if xr.open_dataset(f).sizes.get('iDEPTH', 0) == N_DEPTH_REF),
    None
)
if _ref_file is None:
    raise RuntimeError(f"No post-NCEI file found with iDEPTH={N_DEPTH_REF}")
_ds_ref = xr.open_dataset(_ref_file)
DEPTH_REF = _ds_ref['depth'].values.astype(float)
_ds_ref.close()
print(f"  Reference depth grid: {N_DEPTH_REF} levels, "
      f"{DEPTH_REF[0]:.0f}–{DEPTH_REF[-1]:.0f} m")

files = sorted(PFL_DIR.glob('*.nc'))
print(f"\nFound {len(files)} post-NCEI PFL files")

bin_profiles = defaultdict(list)

for fi, fpath in enumerate(files):
    stem_parts = fpath.stem.split('_')
    prof_type  = stem_parts[4][0] if len(stem_parts) > 4 else '?'
    print(f"[{fi+1}/{len(files)}] {fpath.name}", flush=True)

    try:
        ds = xr.open_dataset(fpath)
    except Exception as e:
        print(f"  skipping: {e}")
        continue

    n_prof = ds.sizes.get('iPROF', 0)
    if n_prof == 0:
        ds.close()
        continue

    yyyymmdd = ds['prof_YYYYMMDD'].values.astype(int)
    years    = yyyymmdd // 10000
    months   = (yyyymmdd % 10000) // 100
    lons     = ds['prof_lon'].values.astype(float)
    lats     = ds['prof_lat'].values.astype(float)
    T_all    = ds['prof_T'].values.astype(float)       # (iPROF, iDEPTH_file)
    S_all    = ds['prof_S'].values.astype(float)
    Tc_all   = ds['prof_Tclim'].values.astype(float)
    Sc_all   = ds['prof_Sclim'].values.astype(float)
    bin_ids  = ds[BIN_ID_VAR].values.astype(float)

    for i in range(n_prof):
        bid = bin_ids[i]
        if np.isnan(bid) or bid <= 0:
            continue
        bin_profiles[int(bid)].append({
            'year':     int(years[i]),
            'month':    int(months[i]),
            'file':     fpath.name,
            'prof_idx': i,
            'date':     int(yyyymmdd[i]),
            'lon':      float(lons[i]),
            'lat':      float(lats[i]),
            'type':     prof_type,
            'T':        pad_to_ref(T_all[i]),
            'S':        pad_to_ref(S_all[i]),
            'Tclim':    pad_to_ref(Tc_all[i]),
            'Sclim':    pad_to_ref(Sc_all[i]),
        })

    ds.close()
    print(f"  {n_prof:,} profiles → {len(bin_profiles):,} bins so far")

print(f"\nAggregating {len(bin_profiles):,} bins...")

# Total calendar months spanned (probability denominator)
all_ym = set()
for profs in bin_profiles.values():
    for p in profs:
        all_ym.add((p['year'], p['month']))
all_ym = sorted(all_ym)
if all_ym:
    y0, m0 = all_ym[0]
    y1, m1 = all_ym[-1]
    total_months = (y1 - y0) * 12 + (m1 - m0) + 1
else:
    total_months = 1

bins_out = {}

for bid, profs in bin_profiles.items():
    ym_groups = defaultdict(list)
    for p in profs:
        ym_groups[(p['year'], p['month'])].append(p)

    months_dict    = {}
    all_T_monthly  = []
    all_S_monthly  = []

    for ym, group in sorted(ym_groups.items()):
        T_mean     = np.nanmean(np.stack([p['T']     for p in group], axis=0), axis=0)
        S_mean     = np.nanmean(np.stack([p['S']     for p in group], axis=0), axis=0)
        Tc_mean    = np.nanmean(np.stack([p['Tclim'] for p in group], axis=0), axis=0)
        Sc_mean    = np.nanmean(np.stack([p['Sclim'] for p in group], axis=0), axis=0)
        all_T_monthly.append(T_mean)
        all_S_monthly.append(S_mean)

        months_dict[ym] = {
            'T':      T_mean,
            'S':      S_mean,
            'Tclim':  Tc_mean,
            'Sclim':  Sc_mean,
            'sources': [{
                'file':     p['file'],
                'prof_idx': p['prof_idx'],
                'date':     p['date'],
                'lon':      p['lon'],
                'lat':      p['lat'],
                'type':     p['type'],
            } for p in group],
        }

    T_long_mean = np.nanmean(np.stack(all_T_monthly, axis=0), axis=0) if all_T_monthly else np.full(N_DEPTH_REF, np.nan)
    S_long_mean = np.nanmean(np.stack(all_S_monthly, axis=0), axis=0) if all_S_monthly else np.full(N_DEPTH_REF, np.nan)

    # Raw ingredients for on-the-fly probability under different denominators.
    # hit_counts[N] = number of months in which this bin received >= N profiles.
    hit_counts = {1: 0, 2: 0, 3: 0}
    for ym, group in ym_groups.items():
        n = len(group)
        if n >= 1: hit_counts[1] += 1
        if n >= 2: hit_counts[2] += 1
        if n >= 3: hit_counts[3] += 1

    bin_yms   = sorted(ym_groups.keys())
    first_ym  = bin_yms[0]
    last_ym   = bin_yms[-1]

    bins_out[bid] = {
        'bin_id':     bid,
        'months':     months_dict,
        'T_mean':     T_long_mean,
        'S_mean':     S_long_mean,
        'hit_counts': hit_counts,
        'first_ym':   first_ym,
        'last_ym':    last_ym,
        # Precomputed day-0 probability retained for backward compatibility
        'prob_1plus': hit_counts[1] / total_months,
    }

result = {
    'depth':             DEPTH_REF,
    'bins':              bins_out,
    'bin_centres':       bin_centres,
    'dataset_first_ym':  (y0, m0) if all_ym else None,
    'dataset_last_ym':   (y1, m1) if all_ym else None,
    'total_months':      total_months,
}

print(f"\nSaving to {OUTPUT_PICKLE}...")
OUTPUT_PICKLE.parent.mkdir(parents=True, exist_ok=True)
with open(OUTPUT_PICKLE, 'wb') as f:
    pickle.dump(result, f, protocol=4)

size_mb = OUTPUT_PICKLE.stat().st_size / 1e6
print(f"Done. {len(bins_out):,} bins saved. File size: {size_mb:.1f} MB")
