"""
``ecco-insitu run``: the whole NCEI chain on a folder of input files.

Steps 01-10 (:mod:`ecco_insitu.pipeline`), then, on the whole run folder, step 11 (spatial
scaling), optionally the step 11 figures and step 12 (compaction), and step 13 (packing into
model-ready files).

- **Settings:** ``ecco_insitu/defaults.yaml``, ``config/local.yaml`` (this machine's data
  locations), ``--config FILE`` and ``--set key=value``, merged in that order.
- **Input and output:** the input folder is ``run.input_dir`` unless ``-i`` is given. The output
  goes to ``run.output_base`` + ``_<YYYYMMDD_HHMMSS>``, so no run overwrites another, unless
  ``-d`` names an exact folder, which must not exist yet.
- **Output layout** (inside the run folder):

  ``yearly/``
      one file per input file (steps 01-10, then 11 and 12)
  ``ecco_arrays/``
      the model-ready profile files (step 13)
  ``documentation/``
      ``run_config.yaml`` (the settings used), ``run_summary.txt`` (what happened to each input
      file), and step 13's manifest and ``PROFILES_SIZE.h`` / ``data.profiles`` values

- **Log:** everything printed, by every step, is also written to
  ``<log_dir>/ncei_log_<YYYYMMDD_HHMMSS>.log`` (``run.log_dir``; default ``z_logs/`` in the
  repository).

Exit codes:

- 0 = all good;
- 2 = bad options or settings (nothing run);
- 3 = finished, but some input files failed (listed in ``run_summary.txt``);
- 1 = a whole-run step (11-13) failed.
"""
import argparse
import contextlib
import datetime
import os
import re
import sys
import threading
from pathlib import Path

from ecco_insitu import config

#: Exit codes of :func:`run_chain` (see the module description).
EXIT_OK, EXIT_STEP_FAILED, EXIT_BAD_SETTINGS, EXIT_FILES_FAILED = 0, 1, 2, 3

_RULE = '=' * 63


@contextlib.contextmanager
def _tee_output(log_path):
    """Send everything written to stdout and stderr, at the OS level (so including C libraries and
    the worker processes, which inherit the file descriptors), to both the console and
    ``log_path``, as ``... 2>&1 | tee log`` does in a shell.

    The copying thread stops when it reads an end marker written after the run, not when the pipe
    closes: helper processes (e.g. multiprocessing's resource tracker) can inherit the pipe and
    keep it open after the run.
    """
    end = b'\0ECCO_INSITU_END_OF_RUN_OUTPUT\0'
    sys.stdout.flush()
    sys.stderr.flush()
    console = os.dup(1)                          # the real console, for the copy thread
    saved_out, saved_err = os.dup(1), os.dup(2)
    r, w = os.pipe()
    os.set_inheritable(r, False)
    os.dup2(w, 1)
    os.dup2(w, 2)

    def copy():
        buf = b''
        with open(log_path, 'wb') as logf:
            while True:
                chunk = os.read(r, 65536)
                if not chunk:
                    break
                buf += chunk
                i = buf.find(end)
                out, buf = (buf[:i], None) if i >= 0 else (buf[:-len(end)], buf[-len(end):])
                os.write(console, out)
                logf.write(out)
                logf.flush()
                if buf is None:
                    break
            if buf:
                os.write(console, buf)
                logf.write(buf)

    t = threading.Thread(target=copy, daemon=True)
    t.start()
    try:
        yield
    finally:
        sys.stdout.flush()
        sys.stderr.flush()
        os.write(w, end)
        os.dup2(saved_out, 1)
        os.dup2(saved_err, 2)
        os.close(saved_out)
        os.close(saved_err)
        t.join(timeout=60)
        os.close(w)
        os.close(r)
        os.close(console)


def _banner(*lines):
    print('\n' + _RULE)
    for line in lines:
        print(line)
    print(_RULE, flush=True)


def default_log_dir(cfg):
    """The folder for run logs: ``run.log_dir`` if set, else ``z_logs/`` in the repository checkout (or
    in the current folder, for a package installed without the repository).

    Parameters
    ----------
    cfg : dict
        Settings from :func:`ecco_insitu.config.load`.

    Returns
    -------
    Path
    """
    if cfg['run'].get('log_dir'):
        return Path(cfg['run']['log_dir'])
    root = config.CONFIG_DIR.parent
    return (root if (root / 'pyproject.toml').exists() else Path.cwd()) / 'z_logs'


def _compact_yearly(yearly_dir):
    """Step 12 on every file under ``yearly_dir``, in place: each file is rewritten compacted and
    renamed ``..._ncei_step_12.nc``. Returns the number of files that failed (left as they were)."""
    from ecco_insitu.steps import step12_profiles_compact as step12
    failed = 0
    for f in sorted(Path(yearly_dir).rglob('*.nc')):
        tmp = f.with_name(f.stem + '_compact.nc')
        try:
            rc = step12.cli([str(f), '-o', str(tmp), '--complevel', '2', '--overwrite'])
        except SystemExit as e:                       # argparse / explicit exits inside the step
            rc = e.code
        if rc:
            print(f'WARNING: step12 failed for {f}; leaving original in place')
            tmp.unlink(missing_ok=True)
            failed += 1
            continue
        new = f.with_name(re.sub(r'__ncei_step_\d+', '__ncei_step_12', f.name))
        tmp.replace(new)
        if new != f:
            f.unlink()
    return failed


def run_chain(input_dir, dest_dir, cfg, n_workers=None, scaling=None, window=None, plots=None,
              compact=None, pack=None):
    """Run the chain; see the module description. Returns an exit code.

    Parameters
    ----------
    input_dir, dest_dir : str or Path
        Input folder; run folder (must not exist).
    cfg : dict
        Settings from :func:`ecco_insitu.config.load`.
    n_workers : int, optional
        Parallel workers (default ``run.n_workers``).
    scaling, plots, compact, pack : bool, optional
        Run step 11, the step 11 figures, step 12, step 13 (defaults: ``run.spatial_scaling``,
        ``run.scaling_plots``, ``run.compact``, ``run.pack``).
    window : str, optional
        Step 11 window, ``monthly`` or ``10day`` (default ``spatial_scaling.window``).
    """
    from ecco_insitu import pipeline
    from ecco_insitu.steps import step11_spatial_scaling as step11
    from ecco_insitu.steps import step13_pack_ecco_arrays as step13

    r = cfg['run']
    scaling = r['spatial_scaling'] if scaling is None else scaling
    plots = r['scaling_plots'] if plots is None else plots
    compact = r['compact'] if compact is None else compact
    pack = r['pack'] if pack is None else pack
    window = window or cfg['spatial_scaling']['window']
    n_workers = r['n_workers'] if n_workers is None else n_workers
    dest_dir = Path(dest_dir)
    yearly, arrays = dest_dir / 'yearly', dest_dir / 'ecco_arrays'

    _banner(f'NCEI chain run started: {datetime.datetime.now():%c}', f'Input dir:  {input_dir}',
            f'Output dir: {dest_dir}')
    try:
        n_failed = pipeline.main(str(yearly), str(input_dir), n_workers=n_workers, cfg=cfg)
    except config.ConfigError as e:
        print(f'config error: {e}', file=sys.stderr)
        return EXIT_BAD_SETTINGS
    if n_failed:
        print(f'WARNING: some input files failed in steps 01-10; see {dest_dir}/documentation/run_summary.txt')
    _banner(f'Steps 01-10 finished: {datetime.datetime.now():%c}', f'Output dir: {dest_dir}')

    status = EXIT_OK
    try:
        if scaling:
            _banner(f'Step 11 (spatial scaling): window={window}')
            step11.main(str(yearly), window_mode=window, n_workers=n_workers)
            _banner(f'Step 11 (spatial scaling) finished: {datetime.datetime.now():%c}')
        if plots:
            from ecco_insitu.diagnostics import report_spatial_scaling
            _banner('Spatial scaling diagnostic plots')
            report_spatial_scaling.main(str(yearly), window_mode=window)
        if compact:
            _banner('Step 12 (compact): rewriting output files in place')
            if _compact_yearly(yearly):
                status = EXIT_STEP_FAILED
            _banner(f'Step 12 (compact) finished: {datetime.datetime.now():%c}')
        if pack and status == EXIT_OK:
            _banner(f'Step 13 (ECCO-ready files): {yearly} -> {arrays}')
            step13.main(str(yearly), str(arrays), n_workers=n_workers, cfg=cfg)
            _banner(f'Step 13 finished: {datetime.datetime.now():%c}')
    except Exception as e:                            # a whole-run step crashed: say which, keep the log
        import traceback
        traceback.print_exc()
        print(f'ERROR: a whole-run step failed: {type(e).__name__}: {e}')
        return EXIT_STEP_FAILED
    if status == EXIT_OK and n_failed:
        print(f'Run finished, but some input files failed: see {dest_dir}/documentation/run_summary.txt')
        return EXIT_FILES_FAILED
    return status


def cli(argv=None):
    """Command-line entry point (``ecco-insitu run``).

    Checks the options and settings, picks the input and output folders, and runs
    :func:`run_chain` with all output also written to the run log.

    Parameters
    ----------
    argv : list of str, optional
        Arguments (default: ``sys.argv[1:]``).

    Returns
    -------
    int
        Exit code; see the module description.
    """
    ap = argparse.ArgumentParser(prog='ecco-insitu run', description=__doc__.split('\n\n')[0],
                                 formatter_class=argparse.RawDescriptionHelpFormatter,
                                 epilog='See the ecco_insitu.run module description for the output layout.')
    ap.add_argument('-i', '--input_dir', help='input folder (default: run.input_dir)')
    ap.add_argument('-d', '--dest_dir',
                    help='exact output folder, must not exist (default: run.output_base_<timestamp>)')
    ap.add_argument('-n', '--n_workers', type=int, help='parallel workers (default: run.n_workers)')
    ap.add_argument('-w', '--window', choices=['monthly', '10day'], help='step 11 windows (spatial_scaling.window)')
    onoff = ap.add_argument_group('steps after 01-10 (defaults from the run.* settings)')
    onoff.add_argument('--no-scaling', dest='scaling', action='store_false', default=None, help='skip step 11')
    onoff.add_argument('--plots', dest='plots', action='store_true', default=None, help='make the step 11 figures')
    onoff.add_argument('--compact', dest='compact', action='store_true', default=None, help='run step 12')
    onoff.add_argument('--no-pack', dest='pack', action='store_false', default=None, help='skip step 13')
    ap.add_argument('--log_dir', help='folder for the run log (default: run.log_dir, else z_logs/)')
    config.add_arguments(ap)
    a = ap.parse_args(argv)
    try:
        cfg = config.from_args(a)
    except config.ConfigError as e:
        print(f'config error: {e}', file=sys.stderr)
        return EXIT_BAD_SETTINGS

    stamp = datetime.datetime.now().strftime('%Y%m%d_%H%M%S')
    input_dir = a.input_dir or cfg['run']['input_dir']
    if not input_dir:
        print('ERROR: no input folder: give -i, or set run.input_dir in config/local.yaml', file=sys.stderr)
        return EXIT_BAD_SETTINGS
    if not Path(input_dir).is_dir():
        print(f'ERROR: input folder does not exist: {input_dir}', file=sys.stderr)
        return EXIT_BAD_SETTINGS
    if a.dest_dir:
        dest = Path(a.dest_dir)
    elif cfg['run']['output_base']:
        dest = Path(f"{cfg['run']['output_base']}_{stamp}")
    else:
        print('ERROR: no output folder: give -d, or set run.output_base in config/local.yaml', file=sys.stderr)
        return EXIT_BAD_SETTINGS
    if dest.exists():
        print(f'ERROR: output folder already exists, refusing to overwrite: {dest}\n'
              '       (pass a different -d, or remove it first)', file=sys.stderr)
        return EXIT_BAD_SETTINGS

    log_dir = Path(a.log_dir) if a.log_dir else default_log_dir(cfg)
    log_dir.mkdir(parents=True, exist_ok=True)
    log_path = log_dir / f'ncei_log_{stamp}.log'
    with _tee_output(log_path):
        print(f'Log file:   {log_path}')
        code = run_chain(input_dir, dest, cfg, n_workers=a.n_workers, scaling=a.scaling, window=a.window,
                         plots=a.plots, compact=a.compact, pack=a.pack)
        _banner(f'NCEI chain run finished: {datetime.datetime.now():%c} (exit code {code})',
                f'Output dir: {dest}', f'Log file:   {log_path}')
    return code


if __name__ == '__main__':
    sys.exit(cli())
