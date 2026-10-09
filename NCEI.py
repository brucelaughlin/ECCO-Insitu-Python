import pdb
import xarray as xr
import numpy as np
import cartopy
import matplotlib.pyplot as plt
import importlib
import sys
import os
import logging
import glob
import io
import contextlib
import tempfile
import concurrent.futures
from pathlib import Path
import argparse
from functools import partial

# Add the directory containing the package to the search path
sys.path.append(os.path.abspath("/Users/brucel/ecco/yip/ECCO-Insitu-Python"))

import tools
import step01
import step02
import step03
import step04
import step05
import step06
import step07
import step08
import step09
import step10


# ==============================================================================
# Worker-process globals: loaded once per worker process at startup.
# On macOS the multiprocessing start method is "spawn" (fresh interpreter),
# so these are not inherited from the main process — the initializer runs in
# every worker and loads the arrays there.
# ==============================================================================
_worker_prebaked_clim_arrays = None


def _worker_init(prebaked_clim_paths):
    """Load pre-baked climatology arrays into a module global, once per worker.

    Keyed by the file's own obs_depth grid, as a tuple of Python floats — the
    same key step03 builds from each input file's prof_depth.
    """
    global _worker_prebaked_clim_arrays
    if prebaked_clim_paths:
        _worker_prebaked_clim_arrays = {}
        for pb_path in prebaked_clim_paths:
            pb_ds = xr.open_dataset(pb_path)
            depth_key = tuple(pb_ds['obs_depth'].values.astype(float).tolist())
            if depth_key in _worker_prebaked_clim_arrays:
                raise ValueError(
                    f"Pre-baked climatologies {Path(pb_path).name} and "
                    f"{_worker_prebaked_clim_arrays[depth_key]['_source_name']} "
                    f"share the same depth grid; remove one.")
            _worker_prebaked_clim_arrays[depth_key] = {
                'prof_T':       pb_ds['potential_T_monthly'].values,
                'prof_S':       pb_ds['S_monthly'].values,
                'lon':          pb_ds['lon'].values,
                'lat':          pb_ds['lat'].values,
                'depths':       pb_ds['obs_depth'].values,
                '_source_name': Path(pb_path).name,
            }
            pb_ds.close()
    else:
        _worker_prebaked_clim_arrays = None


def _normalize_date_time(yyyymmdd, hhmmss):
    """Roll out-of-range YYYYMMDD / HHMMSS values forward into valid ones.

    WOD decimal-time rounding produces SS == 60, MM == 60, HH == 24 (also the
    24:00:00 end-of-day midnight convention) and days past month end (Sep 31).
    Each is treated as an overflow and carried, exactly as prof_date's
    day-number + seconds/86400 arithmetic does:
        042760  -> 042800            235960 day X -> 000000 day X+1
        240000 day X -> 000000 day X+1      Sep 31 -> Oct 1
    Values outside the plausible overflow range (month not 1-12, day not 1-31,
    HH > 24, MM/SS > 60, fill values) are returned unchanged.
    """
    d = np.asarray(yyyymmdd, dtype=np.int64)
    t = np.asarray(hhmmss,   dtype=np.int64)
    y, m, dd = d // 10000, (d // 100) % 100, d % 100
    hh, mi, ss = t // 10000, (t // 100) % 100, t % 100
    ok = ((m >= 1) & (m <= 12) & (dd >= 1) & (dd <= 31) &
          (t >= 0) & (hh <= 24) & (mi <= 60) & (ss <= 60))

    out_d, out_t = d.copy(), t.copy()
    if ok.any():
        month0 = ((y[ok] - 1970).astype('M8[Y]').astype('M8[M]') + (m[ok] - 1))
        ts = (month0.astype('M8[s]')
              + ((dd[ok] - 1) * 86400 + hh[ok] * 3600 + mi[ok] * 60 + ss[ok]))
        day = ts.astype('M8[D]')
        sec = (ts - day.astype('M8[s]')).astype(np.int64)
        yy  = day.astype('M8[Y]').astype(np.int64) + 1970
        mm_ = day.astype('M8[M]').astype(np.int64) % 12 + 1
        dd_ = (day - day.astype('M8[M]').astype('M8[D]')).astype(np.int64) + 1
        out_d[ok] = yy * 10000 + mm_ * 100 + dd_
        out_t[ok] = sec // 3600 * 10000 + (sec // 60) % 60 * 100 + sec % 60
    return out_d, out_t


def _redirect_fd(target_fd):
    """Context manager: redirect stdout at the OS fd level to target_fd.

    contextlib.redirect_stdout only intercepts Python-level sys.stdout.  C
    extensions (netCDF4, HDF5) write directly to fd 1, bypassing it.  This
    context manager duplicates the real fd 1 for later restoration, then
    replaces fd 1 with target_fd for the duration of the block, so ALL output
    — Python and C — goes to the file.
    """
    import contextlib
    saved_fd = os.dup(1)           # save real stdout fd
    os.dup2(target_fd, 1)          # point fd 1 at our temp file
    sys.stdout.flush()
    try:
        yield
    finally:
        sys.stdout.flush()
        os.dup2(saved_fd, 1)       # restore real stdout fd
        os.close(saved_fd)


_redirect_fd = contextlib.contextmanager(_redirect_fd)


def _process_one_file(args):
    """Process a single input file through the full NCEI chain.

    Captures ALL stdout (Python and C-level) to a temp file by redirecting
    fd 1.  Returns the temp file path; the main process reads and prints it
    atomically so the log stays clean even with multiple workers running.

    Returns (file_index, n_total, tmp_log_path, bad_flag).
    """
    (file_dex, n_total, original_file, dest_dir, input_dir,
     ncei_function_kwargs, profile_var_key_set) = args

    original_file = Path(original_file)

    tmp = tempfile.NamedTemporaryFile(mode='w', suffix='.log', delete=False)
    tmp_path = tmp.name
    tmp.close()

    bad_flag = False

    with open(tmp_path, 'w') as log_fh:
        with _redirect_fd(log_fh.fileno()):

            width = len(str(n_total))
            print(f"ncei processing for file {file_dex+1:0{width}}/{n_total}: {original_file}")

            try:
                MITprof_ds = xr.open_dataset(original_file)
                MITprof_ds = MITprof_ds.assign_coords(
                    {dim: np.arange(MITprof_ds.sizes[dim])
                     for dim in MITprof_ds.dims if dim not in MITprof_ds.coords})
            except Exception as e:
                print(f"\n-------- FILE FAILURE --------")
                print(f"Could not open file: {original_file.name}")
                print(f"Exception: {e}")
                print(f"NO OUTPUT FILE WRITTEN for: {original_file.name}")
                print("-----------------------------\n")
                return file_dex, n_total, tmp_path, True

            # Normalize out-of-range dates/times (WOD rounding artefacts, the
            # 24:00:00 midnight convention) by rolling forward, so YYYYMMDD/HHMMSS
            # agree with prof_date.  See _normalize_date_time.
            if 'prof_YYYYMMDD' in MITprof_ds and 'prof_HHMMSS' in MITprof_ds:
                raw_d = MITprof_ds['prof_YYYYMMDD'].values
                raw_t = MITprof_ds['prof_HHMMSS'].values
                fixed_d, fixed_t = _normalize_date_time(raw_d, raw_t)
                n_changed = int(np.count_nonzero((fixed_d != raw_d) | (fixed_t != raw_t)))
                if n_changed:
                    print(f"  normalizing {n_changed} out-of-range prof_YYYYMMDD/prof_HHMMSS values")
                    for name, fixed in (('prof_YYYYMMDD', fixed_d), ('prof_HHMMSS', fixed_t)):
                        MITprof_ds[name] = xr.DataArray(
                            fixed.astype(MITprof_ds[name].dtype),
                            dims=MITprof_ds[name].dims,
                            attrs=MITprof_ds[name].attrs)

            valid_data_dict_list = [tools.collect_valid_data_stats(MITprof_ds, profile_var_key_set)]

            if not MITprof_ds or MITprof_ds.sizes['iPROF'] == 0:
                print("Your profile file may have no valid data; exiting without finishing\n")
                return file_dex, n_total, tmp_path, True

            kw = ncei_function_kwargs
            ncei_function_list = [
                partial(step01.main, grid_dir=kw['grid_dir'], llcN=kw['llcN'], wet_or_all=kw['wet_or_all']),
                partial(step02.main, sphere_bin_dir=kw['sphere_bin_dir'], grid_dir=kw['grid_dir']),
                partial(step03.main, profile_var_key_set=profile_var_key_set,
                        climatology_file=kw['climatology_file'],
                        prebaked_clim_files=_worker_prebaked_clim_arrays),
                partial(step04.main, profile_var_key_set=profile_var_key_set, grid_dir=kw['grid_dir'],
                        sigma_file_dict=kw['sigma_file_dict'],
                        respect_existing_zero_weights=kw['respect_existing_zero_weights'],
                        new_floor_dict=kw['new_floor_dict']),
                partial(step05.main, profile_var_key_set=profile_var_key_set, grid_dir=kw['grid_dir'],
                        apply_gamma_factor=kw['apply_gamma_factor'], llcN=kw['llcN']),
                partial(step06.main, replace_missing_S_with_clim_S=kw['replace_missing_S_with_clim_S']),
                partial(step07.main, profile_var_key_set=profile_var_key_set,
                        exclude_high_latitude_profiles_from_clim_cost=kw['exclude_high_latitude_profiles_from_clim_cost'],
                        dubious_clim_lat_threshold=kw['dubious_clim_lat_threshold']),
                partial(step08.main, profile_var_key_set=profile_var_key_set),
                partial(step09.main, profile_var_key_set=profile_var_key_set),
                partial(step10.main, profile_var_key_set=profile_var_key_set,
                        distance_tolerance=kw['distance_tolerance'],
                        closest_time=kw['closest_time'], method=kw['method']),
            ]

            step_counter = 0
            for ii, fn in enumerate(ncei_function_list):
                step_counter += 1
                print(f"\nstep: {step_counter}")
                try:
                    MITprof_ds = fn(MITprof_ds)
                except Exception as e:
                    print(f"\n-------- FILE FAILURE --------")
                    print(f"NCEI chain CRASHED at step: {step_counter:02d} (of {len(ncei_function_list)})")
                    print(f"File: {original_file.name}")
                    print(f"Exception:\n{e}")
                    print(f"NO OUTPUT FILE WRITTEN for: {original_file.name}")
                    print("continuing to next file")
                    print("-----------------------------\n")
                    bad_flag = True
                    break

                if not MITprof_ds or MITprof_ds.sizes['iPROF'] == 0:
                    print(f"-------- FILE FAILURE --------")
                    print(f"NCEI chain emptied all valid data at step: {step_counter:02d} (of {len(ncei_function_list)})")
                    print(f"File: {original_file.name}")
                    print(f"NO OUTPUT FILE WRITTEN for: {original_file.name}")
                    print("continuing to next file")
                    print("-----------------------------\n")
                    bad_flag = True
                    break

                valid_data_dict_list.append(tools.collect_valid_data_stats(MITprof_ds, profile_var_key_set))
                for prof_key in profile_var_key_set:
                    if prof_key in MITprof_ds:
                        if valid_data_dict_list[0][prof_key]['valid_profile_count'] > 0:
                            print(f"valid {prof_key} profile count / original valid {prof_key} profile count: "
                                   f"{valid_data_dict_list[step_counter][prof_key]['valid_profile_count']}/"
                                   f"{valid_data_dict_list[0][prof_key]['valid_profile_count']} = "
                                   f"{valid_data_dict_list[step_counter][prof_key]['valid_profile_count'] / valid_data_dict_list[0][prof_key]['valid_profile_count']*100:.2f}%")
                for prof_key in profile_var_key_set:
                    if prof_key in MITprof_ds:
                        if valid_data_dict_list[0][prof_key]['valid_data_count'] > 0:
                            print(f"valid {prof_key} data count / original valid {prof_key} data count: "
                                   f"{valid_data_dict_list[step_counter][prof_key]['valid_data_count']}/"
                                   f"{valid_data_dict_list[0][prof_key]['valid_data_count']} = "
                                   f"{valid_data_dict_list[step_counter][prof_key]['valid_data_count'] / valid_data_dict_list[0][prof_key]['valid_data_count']*100:.2f}%")

            if bad_flag or not MITprof_ds or MITprof_ds.sizes['iPROF'] == 0:
                if not bad_flag:
                    print(f"NO OUTPUT FILE WRITTEN for: {original_file.name} (no valid data after all steps)\n")
                return file_dex, n_total, tmp_path, True

            if tools.count_total_survivors_TS(MITprof_ds, profile_var_key_set) > 0:
                for prof_key in profile_var_key_set:
                    if prof_key in MITprof_ds and f'{prof_key}clim' in MITprof_ds and f'{prof_key}weight' in MITprof_ds:
                        MITprof_ds[f'{prof_key}cost'] = (MITprof_ds[prof_key] - MITprof_ds[f'{prof_key}clim'])**2 * MITprof_ds[f'{prof_key}weight']

                if 'prof_lon' in MITprof_ds:
                    lon = MITprof_ds['prof_lon']
                    MITprof_ds['prof_lon'] = xr.where(lon <= 360, ((lon + 180) % 360) - 180, lon)

                print()
                tools.MITprof_write_to_nc(dest_dir, MITprof_ds, len(ncei_function_list), original_file, input_dir)
            else:
                print(f"NO OUTPUT FILE WRITTEN for: {original_file.name} (0 surviving T/S profiles after all steps)\n\n")

    return file_dex, n_total, tmp_path, bad_flag


def NCEI_pipeline(dest_dir, input_dir, n_workers=None):

    input_profile_files = list(Path(input_dir).rglob("*.nc"))
    input_profile_files.sort()

    # ==========================================================================================
    # ========================== START OF NEED PATHS/ PARAMETERS ==============================
    # ==========================================================================================

    profile_var_key_set = {'prof_T', 'prof_S'}

    grid_dir = '/Users/brucel/ecco/yip/sample_data/ecco-insitu/sweet_gdrive/grid_llc90'
    sphere_bin_dir = '/Users/brucel/ecco/yip/sample_data/ecco-insitu/sweet_gdrive/grid_llc90/sphere_point_distribution'
    climatology_file = "/Users/brucel/ecco/yip/woa23_climatology/woa23_decav91C0_TS_clim_potential_T_1deg_fulldepth.nc"

    # Pre-baked climatologies: depth-interpolation done offline per grid, so step03
    # only needs time-blend + lat/lon lookup at runtime.  Every *_prebaked.nc in
    # this directory is used, keyed by the obs_depth grid stored in the file; any
    # input grid without a match falls back to full-depth interpolation (with a
    # log line saying so).  Regenerate with:  python prebake_woa23_climatology.py
    prebaked_clim_dir = '/Users/brucel/ecco/yip/woa23_climatology'

    sigma_file_dict = {
            'prof_T': '/Users/brucel/ecco/yip/sample_data/ecco-insitu/sweet_gdrive/CTD_sigma_TS/Theta_sigma_smoothed_method_02_masked_merged_capped_extrapolated.bin',
            'prof_S': '/Users/brucel/ecco/yip/sample_data/ecco-insitu/sweet_gdrive/CTD_sigma_TS/Salt_sigma_smoothed_method_02_masked_merged_capped_extrapolated.bin',
    }

    llcN = 90
    wet_or_all = 1
    respect_existing_zero_weights = False
    new_floor_dict = {'prof_T': 0, 'prof_S': 0.005}
    apply_gamma_factor = True
    replace_missing_S_with_clim_S = True
    exclude_high_latitude_profiles_from_clim_cost = True
    dubious_clim_lat_threshold = 60
    distance_tolerance = 5e3
    closest_time = 120000
    method = 1

    # ==========================================================================================
    # ========================== END OF NEED PATHS/ PARAMETERS ================================
    # ==========================================================================================

    prebaked_clim_paths = sorted(Path(prebaked_clim_dir).glob('*_prebaked.nc'))
    print(f"[NCEI] {len(prebaked_clim_paths)} pre-baked climatology file(s) in {prebaked_clim_dir}")
    for _pb in prebaked_clim_paths:
        print(f"         {_pb.name}")

    ncei_function_kwargs = dict(
        grid_dir=grid_dir,
        sphere_bin_dir=sphere_bin_dir,
        climatology_file=climatology_file,
        sigma_file_dict=sigma_file_dict,
        llcN=llcN,
        wet_or_all=wet_or_all,
        respect_existing_zero_weights=respect_existing_zero_weights,
        new_floor_dict=new_floor_dict,
        apply_gamma_factor=apply_gamma_factor,
        replace_missing_S_with_clim_S=replace_missing_S_with_clim_S,
        exclude_high_latitude_profiles_from_clim_cost=exclude_high_latitude_profiles_from_clim_cost,
        dubious_clim_lat_threshold=dubious_clim_lat_threshold,
        distance_tolerance=distance_tolerance,
        closest_time=closest_time,
        method=method,
    )

    n_total = len(input_profile_files)

    if n_workers is None:
        import os as _os
        n_workers = min(8, max(1, (_os.cpu_count() or 1) - 1))

    print()
    print(f"NCEI pipeline: {n_total} files, {n_workers} worker(s)")
    print()

    work_items = [
        (file_dex, n_total, str(f), dest_dir, input_dir, ncei_function_kwargs, profile_var_key_set)
        for file_dex, f in enumerate(input_profile_files)
    ]

    def _flush_log(tmp_path):
        """Read a worker's temp log file, print it, and delete it."""
        try:
            with open(tmp_path) as f:
                print(f.read(), end='')
        finally:
            try:
                os.unlink(tmp_path)
            except OSError:
                pass

    if n_workers == 1:
        # Sequential: run initializer inline so worker globals are populated.
        _worker_init(prebaked_clim_paths)
        for item in work_items:
            _, _, tmp_path, _ = _process_one_file(item)
            _flush_log(tmp_path)
    else:
        with concurrent.futures.ProcessPoolExecutor(
            max_workers=n_workers,
            initializer=_worker_init,
            initargs=(prebaked_clim_paths,),
        ) as executor:
            futures = {executor.submit(_process_one_file, item): item[0] for item in work_items}
            for fut in concurrent.futures.as_completed(futures):
                try:
                    _, _, tmp_path, _ = fut.result()
                except Exception as e:
                    print(f"\n[NCEI] Unexpected worker exception: {e}\n")
                    continue
                _flush_log(tmp_path)


def main(dest_dir, input_dir, n_workers=None):
    NCEI_pipeline(dest_dir, input_dir, n_workers=n_workers)


if __name__ == '__main__':

    parser = argparse.ArgumentParser()

    parser.add_argument("-i", "--input_dir", action="store",
                    help="File path to NETCDF files containing MITprofs info.",
                    dest="input_dir", type=str, required=True)

    parser.add_argument("-d", "--dest_dir", action="store",
                help="File path where you would like to store generated NETCDF file.",
                dest="dest_dir", type=str, required=True)

    parser.add_argument("-n", "--n_workers", action="store",
                help="Number of parallel worker processes (default: min(8, cpu_count-1)).",
                dest="n_workers", type=int, default=None)

    args = parser.parse_args()

    main(args.dest_dir, args.input_dir, n_workers=args.n_workers)
