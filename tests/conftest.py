"""Shared pytest setup: make the ecco_insitu package (src/ layout) importable from this checkout,
also without ``pip install -e .``."""
import os
import sys

SRC = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'src'))
if SRC not in sys.path:
    sys.path.insert(0, SRC)
