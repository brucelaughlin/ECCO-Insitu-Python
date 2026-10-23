"""
Load the NCEI chain settings from YAML files.

The settings are merged from, in order (later ones win):

1. ``config/defaults.yaml`` (in the repository; every setting, with comments)
2. ``config/local.yaml`` (this machine's data locations; not committed)
3. the file named by the environment variable ``ECCO_INSITU_CONFIG``, if set
4. each file given with ``--config`` (any subset of the keys)
5. individual overrides, e.g. ``--set chain.new_floor.prof_S=0.006``

Unknown keys and values of the wrong type are errors, so a typo in a config file
can't be silently ignored. Paths are expanded (``~`` and ``${ENV_VARS}``), and relative
paths are taken from ``data_root``.

Typical use::

    import config_loader
    cfg = config_loader.load(['my_run.yaml'])
    grid_dir = cfg['paths']['grid_dir']

Command line, to see the merged settings::

    python config_loader.py [--config my_run.yaml] [--set key.sub=value]
"""
import argparse
import copy
import os
from pathlib import Path

import yaml

REPO = Path(__file__).resolve().parent
CONFIG_DIR = REPO / 'config'
DEFAULTS_FILE = CONFIG_DIR / 'defaults.yaml'
LOCAL_FILE = CONFIG_DIR / 'local.yaml'
ENV_VAR = 'ECCO_INSITU_CONFIG'

# Keys whose values are paths (dotted). Their entries are expanded and resolved against data_root.
PATH_KEYS = (
    'paths.grid_dir', 'paths.sphere_bin_dir', 'paths.climatology_file', 'paths.prebaked_clim_dir',
    'paths.sigma_files.prof_T', 'paths.sigma_files.prof_S',
    'paths.high_lat_sigma_files.prof_T', 'paths.high_lat_sigma_files.prof_S',
    'paths.geodesic_bin_csv',
    'run.input_dir', 'run.output_base',
    'diagnostics.report_dir', 'diagnostics.cost_report_dir',
)

# Paths that steps 01-10 need; checked before a run (see check_required_paths).
CHAIN_REQUIRED_PATHS = (
    'paths.grid_dir', 'paths.sphere_bin_dir', 'paths.climatology_file',
    'paths.sigma_files.prof_T', 'paths.sigma_files.prof_S',
)


class ConfigError(ValueError):
    """A config file is missing, unreadable, has an unknown key, or has a value of the wrong type."""


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
            raise ConfigError(f'{source}: unknown setting "{key}" (see config/defaults.yaml for the valid ones)')
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
    """Value at a dotted key, e.g. ``get(cfg, 'paths.sigma_files.prof_T')``."""
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
        Every setting from ``config/defaults.yaml``, with overrides applied and paths resolved.
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
    """Raise ConfigError listing every required path that is unset or doesn't exist."""
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
                    'WOA23 one with:  python setup_climatology.py --build')
        raise ConfigError('missing input data; set these in config/local.yaml or a --config file '
                          '(see config/local.example.yaml):\n' + '\n'.join(problems) + hint)


def dump(cfg, path):
    """Write the merged settings to ``path`` as YAML, with a header saying where they came from."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    body = {k: v for k, v in cfg.items() if k != '_sources'}
    header = ('# Settings used for this run (merged from, in order):\n'
              + ''.join(f'#   {s}\n' for s in cfg.get('_sources', [])) + '\n')
    path.write_text(header + yaml.safe_dump(body, sort_keys=False, default_flow_style=False))


def add_arguments(parser):
    """Add the standard ``--config`` and ``--set`` options to an argparse parser."""
    parser.add_argument('--config', action='append', default=[], metavar='FILE',
                        help='YAML settings file (repeatable; later files win). See config/defaults.yaml.')
    parser.add_argument('--set', action='append', default=[], metavar='KEY=VALUE', dest='overrides',
                        help='override one setting, e.g. --set chain.new_floor.prof_S=0.006 (repeatable)')


def from_args(args):
    """Settings from parsed ``--config`` / ``--set`` options (see add_arguments)."""
    return load(args.config, args.overrides)


def _shell_value(v):
    """A setting as text for a shell script: true/false as 1/0, null as an empty string."""
    if v is None:
        return ''
    if isinstance(v, bool):
        return '1' if v else '0'
    return str(v)


if __name__ == '__main__':
    ap = argparse.ArgumentParser(description='Print the merged NCEI chain settings.')
    add_arguments(ap)
    ap.add_argument('--get', nargs='+', metavar='KEY',
                    help='print only these dotted settings, one per line (for shell scripts: '
                         'true/false as 1/0, null as an empty line)')
    a = ap.parse_args()
    try:
        c = from_args(a)
        if a.get:
            print('\n'.join(_shell_value(get(c, k)) for k in a.get))
            raise SystemExit(0)
    except (ConfigError, KeyError) as e:
        raise SystemExit(f'config error: {e}')
    print('# merged from: ' + ', '.join(c['_sources']))
    print(yaml.safe_dump({k: v for k, v in c.items() if k != '_sources'}, sort_keys=False))
