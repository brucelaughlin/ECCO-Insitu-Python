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
import copy
import tempfile
import traceback
import concurrent.futures
from pathlib import Path
import argparse
from functools import partial

# Make the chain's modules importable wherever this is run from
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import config_loader
import netCDF4
import tools
import validate_inputs
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


def check_climatologies(climatology_file, prebaked_clim_dir):
    """Check the climatology and choose the pre-baked files to use; return their paths.

    The full-depth climatology must pass :func:`validate_inputs.check_climatology_file`
    (raises ``config_loader.ConfigError`` otherwise, before any file is processed).
    A pre-baked file is used only if it passes its own check **and** its ``source_file``
    attribute names this climatology; pre-baked files made from another climatology are
    ignored with a message, so changing ``paths.climatology_file`` can't silently mix two
    climatologies.
    """
    rep = validate_inputs.check_climatology_file(climatology_file)
    if not rep.ok:
        raise config_loader.ConfigError(
            f'paths.climatology_file is not usable: {climatology_file}\n  ' + '\n  '.join(rep.lines()) +
            '\nSee docs/input_contract.md, or build one with: python build_woa23_climatology.py')
    for w in rep.warnings:
        print(f"[NCEI] climatology WARNING: {w}")
    print(f"[NCEI] climatology: {climatology_file}")
    found = sorted(Path(prebaked_clim_dir).glob('*_prebaked.nc')) if prebaked_clim_dir else []
    use = []
    for pb in found:
        prep = validate_inputs.check_climatology_file(pb, prebaked=True)
        with netCDF4.Dataset(pb) as d:
            src = d.getncattr('source_file') if 'source_file' in d.ncattrs() else None
        if not prep.ok:
            print(f"[NCEI] pre-baked climatology IGNORED (not usable): {pb.name}: {'; '.join(prep.errors)}")
        elif src is None or Path(src).resolve() != Path(climatology_file).resolve():
            print(f"[NCEI] pre-baked climatology IGNORED (made from a different climatology, {src}): {pb.name}")
        else:
            use.append(pb)
    print(f"[NCEI] {len(use)} pre-baked climatology file(s) in {prebaked_clim_dir}")
    for pb in use:
        print(f"         {pb.name}")
    return use


def _validate_input_files(files, input_dir, variables):
    """Check every input file (validate_inputs.check_profile_file); return (files to process, reports)."""
    reports = [validate_inputs.check_profile_file(f, sorted(variables)) for f in files]
    print()
    print("[NCEI] input file check (see docs/input_contract.md):")
    for line in validate_inputs.summarize(reports, input_dir):
        print(f"  {line}")
    return [Path(r.path) for r in reports if r.ok], reports


def _write_run_summary(run_root, input_dir, input_reports, results):
    """Write <run>/documentation/run_summary.txt: per input file, what happened to it."""
    status = {r[0]: r for r in results}
    lines = ['NCEI chain run summary (steps 01-10)', f'input folder: {input_dir}', '']
    n_skip = sum(not r.ok for r in input_reports)
    n_fail = sum(1 for r in results if r[1] == 'failed')
    n_empty = sum(1 for r in results if r[1] == 'no data left')
    n_ok = sum(1 for r in results if r[1] == 'written')
    lines += [f'files found:                  {len(input_reports)}',
              f'skipped by the input check:   {n_skip}',
              f'written:                      {n_ok}',
              f'no data left after the steps: {n_empty}',
              f'failed:                       {n_fail}', '']
    for rep in input_reports:
        rel = str(Path(rep.path).relative_to(input_dir))
        if not rep.ok:
            lines.append(f'SKIPPED  {rel}')
            lines += [f'           {x}' for x in rep.lines()]
            continue
        _, state, detail = status.get(rep.path, (rep.path, 'not run', ''))
        lines.append(f'{state.upper():8s} {rel}' + (f'  ({detail})' if detail else ''))
        lines += [f'           WARNING {w}' for w in rep.warnings]
    path = Path(run_root) / 'documentation' / 'run_summary.txt'
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text('\n'.join(lines) + '\n')
    print(f"\n[NCEI] {n_ok} written, {n_empty} with no data left, {n_fail} failed, {n_skip} skipped by the input check;"
          f" details in {path}")
    return n_fail


def _process_one_file(args):
    """Process a single input file through the full NCEI chain.

    Captures ALL stdout (Python and C-level) to a temp file by redirecting
    fd 1.  Returns the temp file path; the main process reads and prints it
    atomically so the log stays clean even with multiple workers running.

    Returns (file_index, n_total, tmp_log_path, state, detail): state is 'written',
    'no data left' (QC removed every profile; not an error) or 'failed' (an exception;
    ``detail`` says at which step and why).
    """
    (file_dex, n_total, original_file, dest_dir, input_dir,
     ncei_function_kwargs, profile_var_key_set) = args

    original_file = Path(original_file)

    tmp = tempfile.NamedTemporaryFile(mode='w', suffix='.log', delete=False)
    tmp_path = tmp.name
    tmp.close()

    bad_flag = False
    state, detail = 'failed', ''

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
                return file_dex, n_total, tmp_path, 'failed', f'could not be opened: {e}'

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
                return file_dex, n_total, tmp_path, 'no data left', 'no profiles in the file'

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
                        high_lat_sigma_file_dict=kw['high_lat_sigma_file_dict'],
                        high_lat_sigma_min_lat=kw['high_lat_sigma_min_lat'],
                        high_lat_sigma_max_dist_km=kw['high_lat_sigma_max_dist_km'],
                        southern_sigma_const_dict=kw['southern_sigma_const_dict'],
                        new_floor_dict=kw['new_floor_dict']),
                partial(step05.main, replace_missing_S_with_clim_S=kw['replace_missing_S_with_clim_S'],
                        convert_insitu_T_to_potential=kw['convert_insitu_T_to_potential']),
                partial(step06.main, profile_var_key_set=profile_var_key_set,
                        exclude_high_latitude_profiles_from_clim_cost=kw['exclude_high_latitude_profiles_from_clim_cost'],
                        dubious_clim_lat_threshold=kw['dubious_clim_lat_threshold'],
                        clim_cost_qc_with_high_lat_sigma=kw['clim_cost_qc_with_high_lat_sigma']),
                partial(step07.main, profile_var_key_set=profile_var_key_set),
                partial(step08.main, profile_var_key_set=profile_var_key_set,
                        remove_extraneous_depth_levels=kw['remove_extraneous_depth_levels']),
                partial(step09.main, profile_var_key_set=profile_var_key_set,
                        distance_tolerance=kw['distance_tolerance'],
                        closest_time=kw['closest_time'], method=kw['method']),
                partial(step10.main, profile_var_key_set=profile_var_key_set, grid_dir=kw['grid_dir'],
                        apply_gamma_factor=kw['apply_gamma_factor'], llcN=kw['llcN']),
            ]

            step_counter = 0
            t_conversion = None              # step05's decision, kept for the output file
            for ii, fn in enumerate(ncei_function_list):
                step_counter += 1
                print(f"\nstep: {step_counter}")
                try:
                    MITprof_ds = fn(MITprof_ds)
                    if 'prof_T' in MITprof_ds and 'T_conversion' in MITprof_ds['prof_T'].attrs:
                        t_conversion = MITprof_ds['prof_T'].attrs['T_conversion']
                except Exception as e:
                    print(f"\n-------- FILE FAILURE --------")
                    print(f"NCEI chain CRASHED at step: {step_counter:02d} (of {len(ncei_function_list)})")
                    print(f"File: {original_file.name}")
                    print(f"Exception ({type(e).__name__}):\n{e}")
                    print(traceback.format_exc())
                    print(f"NO OUTPUT FILE WRITTEN for: {original_file.name}")
                    print("continuing to next file")
                    print("-----------------------------\n")
                    bad_flag = True
                    state, detail = 'failed', f'step {step_counter:02d}: {type(e).__name__}: {e}'
                    break

                if not MITprof_ds or MITprof_ds.sizes['iPROF'] == 0:
                    print(f"-------- FILE FAILURE --------")
                    print(f"NCEI chain emptied all valid data at step: {step_counter:02d} (of {len(ncei_function_list)})")
                    print(f"File: {original_file.name}")
                    print(f"NO OUTPUT FILE WRITTEN for: {original_file.name}")
                    print("continuing to next file")
                    print("-----------------------------\n")
                    bad_flag = True
                    state, detail = 'no data left', f'all profiles removed by step {step_counter:02d}'
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
                    state, detail = 'no data left', 'no valid data after all steps'
                return file_dex, n_total, tmp_path, state, detail

            if tools.count_total_survivors_TS(MITprof_ds, profile_var_key_set) > 0:
                for prof_key in profile_var_key_set:
                    if prof_key in MITprof_ds and f'{prof_key}clim' in MITprof_ds and f'{prof_key}weight' in MITprof_ds:
                        MITprof_ds[f'{prof_key}cost'] = (MITprof_ds[prof_key] - MITprof_ds[f'{prof_key}clim'])**2 * MITprof_ds[f'{prof_key}weight']

                if 'prof_lon' in MITprof_ds:
                    lon = MITprof_ds['prof_lon']
                    MITprof_ds['prof_lon'] = xr.where(lon <= 360, ((lon + 180) % 360) - 180, lon)

                # The steps drop variable attributes; say what prof_T is (potential
                # temperature: step05 converts in-situ T, and leaves potential T alone).
                if 'prof_T' in MITprof_ds:
                    MITprof_ds['prof_T'].attrs.update(long_name='potential temperature', units='degree C')
                    if t_conversion:          # how step05 decided (metadata, climatology check, default, forced)
                        MITprof_ds['prof_T'].attrs['T_conversion'] = t_conversion

                # Chain-internal only (steps 01 and 10 use it to look up grid-cell values); not saved.
                MITprof_ds = MITprof_ds.drop_vars('profile_flattened_monotonic_grid_indices', errors='ignore')

                print()
                tools.MITprof_write_to_nc(dest_dir, MITprof_ds, len(ncei_function_list), original_file, input_dir)
                state, detail = 'written', f'{MITprof_ds.sizes["iPROF"]} profiles'
            else:
                print(f"NO OUTPUT FILE WRITTEN for: {original_file.name} (0 surviving T/S profiles after all steps)\n\n")
                state, detail = 'no data left', '0 surviving T/S data after all steps'

    return file_dex, n_total, tmp_path, state, detail


def NCEI_pipeline(dest_dir, input_dir, n_workers=None, cfg=None):
    """Run steps 01-10 on every ``*.nc`` file under ``input_dir``, writing to ``dest_dir``.

    Parameters
    ----------
    dest_dir : str or Path
        Output folder; files keep their sub-folder layout from ``input_dir``.
    input_dir : str or Path
        Folder of input profile files, searched recursively.
    n_workers : int, optional
        Parallel worker processes; default ``run.n_workers`` from the settings, else
        min(8, CPUs - 1).
    cfg : dict, optional
        Settings from :func:`config_loader.load`; default: the standard files
        (``config/defaults.yaml`` + ``config/local.yaml``).

    Returns
    -------
    int
        Number of input files that failed with an error (0 for a clean run). Files skipped by the
        input check and files with no data left after QC are not counted; all are listed in
        ``<run>/documentation/run_summary.txt``.
    """
    if cfg is None:
        cfg = config_loader.load()
    config_loader.check_required_paths(cfg)

    input_profile_files = list(Path(input_dir).rglob("*.nc"))
    input_profile_files.sort()

    paths, ch, hl = cfg['paths'], cfg['chain'], cfg['high_latitude_sigma']
    profile_var_key_set = set(ch['variables'])

    grid_dir = paths['grid_dir']
    sphere_bin_dir = paths['sphere_bin_dir']
    climatology_file = paths['climatology_file']
    # Pre-baked climatologies (optional): every *_prebaked.nc in this folder is used, keyed
    # by its obs_depth grid; an input grid without a match falls back to full-depth
    # interpolation. Regenerate with:  python prebake_woa23_climatology.py
    prebaked_clim_dir = paths['prebaked_clim_dir']
    sigma_file_dict = dict(paths['sigma_files'])

    # High-latitude uncertainties (step 04), used instead of the CTD sigma poleward of
    # min_lat (prof_Vuncertainty_source: 0 = CTD sigma, 1 = Arctic file, 2 = southern constant).
    use_high_lat_sigma = hl['enabled']
    high_lat_sigma_file_dict = (
        {v: (paths['high_lat_sigma_files'][v], hl['north_file_vars'][v]) for v in ('prof_T', 'prof_S')}
        if use_high_lat_sigma else None)
    if use_high_lat_sigma:
        config_loader.check_required_paths(cfg, ('paths.high_lat_sigma_files.prof_T',
                                                 'paths.high_lat_sigma_files.prof_S'))
    high_lat_sigma_min_lat = hl['min_lat']
    high_lat_sigma_max_dist_km = hl['max_dist_km']
    southern_sigma_const_dict = dict(hl['south_constant']) if use_high_lat_sigma and hl['south_constant'] else None

    llcN = ch['llcN']
    wet_or_all = ch['wet_or_all']
    respect_existing_zero_weights = ch['respect_existing_zero_weights']
    new_floor_dict = dict(ch['new_floor'])
    apply_gamma_factor = ch['apply_gamma_factor']
    replace_missing_S_with_clim_S = ch['replace_missing_S_with_clim_S']
    convert_insitu_T_to_potential = ch['convert_insitu_T_to_potential']
    exclude_high_latitude_profiles_from_clim_cost = ch['exclude_high_latitude_profiles_from_clim_cost']
    dubious_clim_lat_threshold = ch['dubious_clim_lat_threshold']
    clim_cost_qc_with_high_lat_sigma = ch['clim_cost_qc_with_high_lat_sigma']
    remove_extraneous_depth_levels = ch['remove_extraneous_depth_levels']
    distance_tolerance = ch['distance_tolerance']
    closest_time = ch['closest_time']
    method = ch['method']
    if n_workers is None:
        n_workers = cfg['run']['n_workers']

    # Record the settings used with the run: <run>/documentation/run_config.yaml
    # (<run> is the parent of dest_dir when dest_dir is <run>/yearly).
    run_root = Path(dest_dir).parent if Path(dest_dir).name == 'yearly' else Path(dest_dir)
    used = copy.deepcopy(cfg)
    used['run']['input_dir'] = str(Path(input_dir).resolve())
    config_loader.dump(used, run_root / 'documentation' / 'run_config.yaml')

    prebaked_clim_paths = check_climatologies(climatology_file, prebaked_clim_dir)

    # Check every input file before processing any (see docs/input_contract.md): files that
    # can't be processed are skipped with the reason; warnings are listed but don't stop a file.
    input_profile_files, input_reports = _validate_input_files(input_profile_files, input_dir, profile_var_key_set)

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
        convert_insitu_T_to_potential=convert_insitu_T_to_potential,
        exclude_high_latitude_profiles_from_clim_cost=exclude_high_latitude_profiles_from_clim_cost,
        dubious_clim_lat_threshold=dubious_clim_lat_threshold,
        clim_cost_qc_with_high_lat_sigma=clim_cost_qc_with_high_lat_sigma,
        high_lat_sigma_file_dict=high_lat_sigma_file_dict,
        high_lat_sigma_min_lat=high_lat_sigma_min_lat,
        high_lat_sigma_max_dist_km=high_lat_sigma_max_dist_km,
        southern_sigma_const_dict=southern_sigma_const_dict,
        remove_extraneous_depth_levels=remove_extraneous_depth_levels,
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
    if not work_items:
        print("[NCEI] no input file passed the input check; nothing to process")

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

    results = []                     # (input path, state, detail) per processed file
    if n_workers == 1:
        # Sequential: run initializer inline so worker globals are populated.
        _worker_init(prebaked_clim_paths)
        for item in work_items:
            _, _, tmp_path, state, detail = _process_one_file(item)
            results.append((item[2], state, detail))
            _flush_log(tmp_path)
    else:
        with concurrent.futures.ProcessPoolExecutor(
            max_workers=n_workers,
            initializer=_worker_init,
            initargs=(prebaked_clim_paths,),
        ) as executor:
            futures = {executor.submit(_process_one_file, item): item[0] for item in work_items}
            for fut in concurrent.futures.as_completed(futures):
                path = work_items[futures[fut]][2]
                try:
                    _, _, tmp_path, state, detail = fut.result()
                except Exception as e:
                    print(f"\n[NCEI] Unexpected worker exception: {e}\n")
                    results.append((path, 'failed', f'worker crashed: {type(e).__name__}: {e}'))
                    continue
                results.append((path, state, detail))
                _flush_log(tmp_path)

    n_failed = _write_run_summary(run_root, input_dir, input_reports, results)
    return n_failed


def main(dest_dir, input_dir, n_workers=None, cfg=None):
    """Run steps 01-10; see :func:`NCEI_pipeline`. Returns the number of files that failed."""
    return NCEI_pipeline(dest_dir, input_dir, n_workers=n_workers, cfg=cfg)


if __name__ == '__main__':

    parser = argparse.ArgumentParser()

    parser.add_argument("-i", "--input_dir", action="store",
                    help="Folder of input profile files (default: run.input_dir in the settings).",
                    dest="input_dir", type=str, default=None)

    parser.add_argument("-d", "--dest_dir", action="store",
                help="File path where you would like to store generated NETCDF file.",
                dest="dest_dir", type=str, required=True)

    parser.add_argument("-n", "--n_workers", action="store",
                help="Number of parallel worker processes (default: min(8, cpu_count-1)).",
                dest="n_workers", type=int, default=None)

    config_loader.add_arguments(parser)

    args = parser.parse_args()
    try:
        cfg = config_loader.from_args(args)
        input_dir = args.input_dir or cfg['run']['input_dir']
        if not input_dir:
            raise config_loader.ConfigError('no input folder: give -i or set run.input_dir')
    except config_loader.ConfigError as e:
        sys.exit(f"config error: {e}")
    try:
        n_failed = main(args.dest_dir, input_dir, n_workers=args.n_workers, cfg=cfg)
    except config_loader.ConfigError as e:
        sys.exit(f"config error: {e}")
    # Exit code 3 = the run finished, but some files failed (see documentation/run_summary.txt).
    sys.exit(3 if n_failed else 0)
