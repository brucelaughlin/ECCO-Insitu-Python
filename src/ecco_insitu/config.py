"""
Load the NCEI chain settings from YAML files.

The settings are merged from, in order (later ones win):

1. ``src/ecco_insitu/defaults.yaml`` (shipped with the package; every setting, with comments)
2. ``config/local.yaml`` in the repository checkout (this machine's data locations; not committed;
   ``~/.config/ecco_insitu/local.yaml`` for a pip-installed package)
3. the file named by the environment variable ``ECCO_INSITU_CONFIG``, if set
4. each file given with ``--config`` (any subset of the keys)
5. individual overrides, e.g. ``--set chain.new_floor.prof_S=0.006``

Unknown keys and values of the wrong type are errors, so a typo in a config file
can't be silently ignored. Paths are expanded (``~`` and ``${ENV_VARS}``), and relative
paths are taken from ``data_root``.

Typical use::

    import ecco_insitu.config as config_loader
    cfg = config_loader.load(['my_run.yaml'])
    grid_dir = cfg['paths']['grid_dir']

Command line, to see the merged settings::

    ecco-insitu config [--config my_run.yaml] [--set key.sub=value]
"""
import argparse
import copy
import os
import sys
from pathlib import Path

import yaml

PACKAGE_DIR = Path(__file__).resolve().parent
#: Every setting with its default; shipped inside the package.
DEFAULTS_FILE = PACKAGE_DIR / 'defaults.yaml'
#: Repository checkout's config folder (holds local.example.yaml, and local.yaml on this machine);
#: the package is at <repo>/src/ecco_insitu.
CONFIG_DIR = PACKAGE_DIR.parent.parent / 'config'
ENV_VAR = 'ECCO_INSITU_CONFIG'


def _local_file():
    """This machine's settings file: ``config/local.yaml`` in the repository checkout if that
    folder exists, else ``~/.config/ecco_insitu/local.yaml`` (for a pip-installed package)."""
    if CONFIG_DIR.is_dir():
        return CONFIG_DIR / 'local.yaml'
    return Path(os.environ.get('XDG_CONFIG_HOME', Path.home() / '.config')) / 'ecco_insitu' / 'local.yaml'


LOCAL_FILE = _local_file()

# Keys whose values are paths (dotted). Their entries are expanded and resolved against data_root.
PATH_KEYS = (
    'paths.grid_dir', 'paths.sphere_bin_dir', 'paths.climatology_file', 'paths.prebaked_clim_dir',
    'paths.sigma_files.prof_T', 'paths.sigma_files.prof_S',
    'paths.high_lat_sigma_files.prof_T', 'paths.high_lat_sigma_files.prof_S',
    'paths.geodesic_bin_csv', 'paths.wod_std_dir',
    'run.input_dir', 'run.output_base', 'run.log_dir',
    'diagnostics.report_dir', 'diagnostics.cost_report_dir',
)

# Paths that steps 01-10 need; checked before a run (see check_required_paths).
CHAIN_REQUIRED_PATHS = (
    'paths.grid_dir', 'paths.sphere_bin_dir', 'paths.climatology_file',
    'paths.sigma_files.prof_T', 'paths.sigma_files.prof_S',
)


class ConfigError(ValueError):
    """A settings problem: a config file is missing or unreadable, has an unknown key or a value of the
    wrong type, or a required data path is unset or doesn't exist.
    """


def _read_yaml(path):
    path = Path(path)
    if not path.is_file():
        raise ConfigError(f'config file not found: {path}')
    try:
        data = yaml.safe_load(path.read_text())
    except yaml.YAMLError as e:
        raise ConfigError(f'{path}: not valid YAML: {e}') from e
    if data is None:
        return {}
    if not isinstance(data, dict):
        raise ConfigError(f'{path}: the top level must be a mapping of settings')
    return data


def _type_ok(default, value):
    """Is ``value`` acceptable where the default is ``default``? null defaults accept anything."""
    if default is None or value is None:
        return True
    if isinstance(default, bool):
        return isinstance(value, bool)
    if isinstance(default, (int, float)) and not isinstance(default, bool):
        return isinstance(value, (int, float)) and not isinstance(value, bool)
    if isinstance(default, str):
        # a string default ('auto', 'monthly', ...) may be replaced by another string, a bool or a number
        return isinstance(value, (str, bool, int, float))
    if isinstance(default, list):
        return isinstance(value, list)
    if isinstance(default, dict):
        return isinstance(value, dict)
    return True


def _merge(base, over, source, prefix=''):
    """Merge mapping ``over`` into ``base`` in place, rejecting keys that aren't in ``base``."""
    for k, v in over.items():
        key = f'{prefix}{k}'
        if k not in base:
            raise ConfigError(f'{source}: unknown setting "{key}" (see src/ecco_insitu/defaults.yaml for the valid ones)')
        if isinstance(base[k], dict) and isinstance(v, dict):
            _merge(base[k], v, source, key + '.')
        elif isinstance(base[k], dict) and v is not None:
            raise ConfigError(f'{source}: "{key}" must be a mapping of settings, not {v!r}')
        elif not _type_ok(base[k], v):
            raise ConfigError(f'{source}: "{key}" = {v!r} has the wrong type '
                              f'(expected {type(base[k]).__name__}, like the default {base[k]!r})')
        else:
            base[k] = v


def _parse_override(text):
    """``a.b.c=value`` -> nested dict, value parsed as YAML (so numbers, true/false, null, lists work)."""
    if '=' not in text:
        raise ConfigError(f'--set {text!r}: expected key.sub=value')
    key, val = text.split('=', 1)
    try:
        value = yaml.safe_load(val)
    except yaml.YAMLError as e:
        raise ConfigError(f'--set {text!r}: value is not valid YAML: {e}') from e
    out = cur = {}
    parts = key.strip().split('.')
    for p in parts[:-1]:
        cur[p] = {}
        cur = cur[p]
    cur[parts[-1]] = value
    return out


def get(cfg, dotted):
    """Value at a dotted key.

    Parameters
    ----------
    cfg : dict
        Settings from :func:`load`.
    dotted : str
        E.g. ``'paths.sigma_files.prof_T'``.

    Returns
    -------
    object

    Raises
    ------
    KeyError
        If the key doesn't exist.
    """
    cur = cfg
    for p in dotted.split('.'):
        cur = cur[p]
    return cur


def _set(cfg, dotted, value):
    parts = dotted.split('.')
    cur = cfg
    for p in parts[:-1]:
        cur = cur[p]
    cur[parts[-1]] = value


def _resolve_paths(cfg):
    root = cfg.get('data_root')
    root = Path(os.path.expandvars(os.path.expanduser(root))) if root else None
    if root is not None:
        cfg['data_root'] = str(root)
    for dotted in PATH_KEYS:
        v = get(cfg, dotted)
        if v is None:
            continue
        p = Path(os.path.expandvars(os.path.expanduser(str(v))))
        if not p.is_absolute():
            if root is None:
                raise ConfigError(f'"{dotted}" = {v!r} is a relative path, but data_root is not set')
            p = root / p
        _set(cfg, dotted, str(p))


def load(config_files=(), overrides=(), use_local=True, use_env=True):
    """Merged settings as a nested dict.

    Parameters
    ----------
    config_files : sequence of str or Path
        Extra config files, applied in order after the defaults, ``config/local.yaml`` and
        ``$ECCO_INSITU_CONFIG``.
    overrides : sequence of str
        ``key.sub=value`` strings, applied last.
    use_local : bool
        Read ``config/local.yaml`` if it exists.
    use_env : bool
        Read the file named by ``$ECCO_INSITU_CONFIG`` if it is set.

    Returns
    -------
    dict
        Every setting from ``src/ecco_insitu/defaults.yaml``, with overrides applied and paths resolved.
        ``cfg['_sources']`` lists the files that were merged, in order.

    Raises
    ------
    ConfigError
        For a missing or unreadable file, an unknown key, or a value of the wrong type.
    """
    cfg = copy.deepcopy(_read_yaml(DEFAULTS_FILE))
    sources = [str(DEFAULTS_FILE)]
    layers = []
    if use_local and LOCAL_FILE.is_file():
        layers.append(LOCAL_FILE)
    if use_env and os.environ.get(ENV_VAR):
        layers.append(Path(os.environ[ENV_VAR]))
    layers.extend(Path(f) for f in config_files)
    for f in layers:
        _merge(cfg, _read_yaml(f), str(f))
        sources.append(str(f))
    for o in overrides:
        _merge(cfg, _parse_override(o), f'--set {o}')
        sources.append(f'--set {o}')
    _resolve_paths(cfg)
    cfg['_sources'] = sources
    return cfg


def check_required_paths(cfg, keys=CHAIN_REQUIRED_PATHS):
    """Check that the given path settings are set and exist.

    Parameters
    ----------
    cfg : dict
        Settings from :func:`load`.
    keys : sequence of str
        Dotted path settings to check (default: those steps 01-10 need, ``CHAIN_REQUIRED_PATHS``).

    Raises
    ------
    ConfigError
        Listing every path that is unset or doesn't exist; for a missing climatology the message also
        says how to build one (``ecco-insitu climatology --build``).
    """
    problems = []
    for k in keys:
        v = get(cfg, k)
        if v is None:
            problems.append(f'  {k}: not set')
        elif not Path(v).exists():
            problems.append(f'  {k}: does not exist: {v}')
    if problems:
        hint = ''
        if any(p.startswith('  paths.climatology_file') for p in problems):
            hint = ('\nFor the climatology: point paths.climatology_file at your own file, or build the standard '
                    'WOA23 one with:  ecco-insitu climatology --build')
        raise ConfigError('missing input data; set these in config/local.yaml or a --config file '
                          '(see config/local.example.yaml):\n' + '\n'.join(problems) + hint)


def dump(cfg, path):
    """Write the merged settings to ``path`` as YAML, with a header listing the files they came from.

    Parameters
    ----------
    cfg : dict
        Settings from :func:`load`.
    path : str or Path
        Output file (its folder is created if needed).
    """
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    body = {k: v for k, v in cfg.items() if k != '_sources'}
    import ecco_insitu   # here, not at the top: the package's __init__ must stay import-free
    header = (f'# Settings used for this run, by {ecco_insitu.provenance()} (merged from, in order):\n'
              + ''.join(f'#   {s}\n' for s in cfg.get('_sources', [])) + '\n')
    path.write_text(header + yaml.safe_dump(body, sort_keys=False, default_flow_style=False))


def add_arguments(parser):
    """Add the standard ``--config FILE`` and ``--set KEY=VALUE`` options (both repeatable) to an
    ``argparse`` parser.

    Parameters
    ----------
    parser : argparse.ArgumentParser
    """
    parser.add_argument('--config', action='append', default=[], metavar='FILE',
                        help='YAML settings file (repeatable; later files win). See src/ecco_insitu/defaults.yaml.')
    parser.add_argument('--set', action='append', default=[], metavar='KEY=VALUE', dest='overrides',
                        help='override one setting, e.g. --set chain.new_floor.prof_S=0.006 (repeatable)')


def from_args(args):
    """Settings from parsed ``--config`` / ``--set`` options (see :func:`add_arguments`).

    Parameters
    ----------
    args : argparse.Namespace

    Returns
    -------
    dict
        As :func:`load`.

    Raises
    ------
    ConfigError
    """
    return load(args.config, args.overrides)


def _shell_value(v):
    """A setting as text for a shell script: true/false as 1/0, null as an empty string."""
    if v is None:
        return ''
    if isinstance(v, bool):
        return '1' if v else '0'
    return str(v)


def cli(argv=None):
    """Command-line entry point (``ecco-insitu config``): print the merged settings, or single values with
    ``--get KEY ...`` (one per line; true/false as 1/0, null as an empty line).

    Parameters
    ----------
    argv : list of str, optional
        Arguments (default: ``sys.argv[1:]``).

    Returns
    -------
    int
        0; exits with a message for a settings error.
    """
    ap = argparse.ArgumentParser(description='Print the merged NCEI chain settings.')
    add_arguments(ap)
    ap.add_argument('--get', nargs='+', metavar='KEY',
                    help='print only these dotted settings, one per line (for shell scripts: '
                         'true/false as 1/0, null as an empty line)')
    a = ap.parse_args(argv)
    try:
        c = from_args(a)
        if a.get:
            print('\n'.join(_shell_value(get(c, k)) for k in a.get))
            return 0
    except (ConfigError, KeyError) as e:
        raise SystemExit(f'config error: {e}')
    print('# merged from: ' + ', '.join(c['_sources']))
    print(yaml.safe_dump({k: v for k, v in c.items() if k != '_sources'}, sort_keys=False))
    return 0


if __name__ == '__main__':
    sys.exit(cli())
