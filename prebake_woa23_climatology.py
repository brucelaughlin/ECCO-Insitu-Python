"""Compatibility wrapper: this script moved to ecco_insitu/climatology/prebake.py.

Same as:  ecco-insitu climatology --prebake   (uses run.input_dir and paths.climatology_file), or,
with explicit --source / --input_dir / --out_dir:  python -m ecco_insitu.climatology.prebake ...
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ecco_insitu.climatology.prebake import cli  # noqa: E402

if __name__ == '__main__':
    sys.exit(cli())
