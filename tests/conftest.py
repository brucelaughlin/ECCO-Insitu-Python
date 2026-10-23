"""Shared pytest setup: make the chain modules importable from the repo root."""
import os
import sys

REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
if REPO not in sys.path:
    sys.path.insert(0, REPO)
