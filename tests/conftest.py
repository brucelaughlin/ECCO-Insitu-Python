"""Shared pytest setup: make the ecco_insitu package (src/ layout) importable from this checkout,
also without ``pip install -e .``."""
import os
import sys

SRC = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'src'))
if SRC not in sys.path:
    sys.path.insert(0, SRC)


def _chain_data_available():
    """True if the settings (config/local.yaml or $ECCO_INSITU_CONFIG) point at the chain's data."""
    try:
        from ecco_insitu import config
        config.check_required_paths(config.load())
        return True
    except Exception:
        return False


#: Whether the full chain can run here (grid, climatology, sigma files found via the settings).
#: Tests that run the chain are skipped when it is False (e.g. on GitHub's servers).
CHAIN_DATA = _chain_data_available()


def climatology_file():
    """The configured climatology file, or None if not set / not present."""
    try:
        from ecco_insitu import config
        p = config.load()['paths']['climatology_file']
        return p if p and os.path.exists(p) else None
    except Exception:
        return None
