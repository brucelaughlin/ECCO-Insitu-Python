"""
Logging for the chain's controller and tools.

Messages go through the standard :mod:`logging` module, under the ``ecco_insitu`` logger, so a
caller can redirect, filter or silence them (e.g. ``logging.getLogger('ecco_insitu').setLevel
(logging.WARNING)``). By default they are printed to stdout as plain text, exactly as before, so
the chain's log files look the same.

The per-file output of steps 01-10 (``step: N``, ``valid prof_T profile count ...``) is
captured per file and printed as one block; diagnostics (``ecco_insitu.diagnostics.costs_scaling``)
parse those lines, so their format is fixed and they stay ordinary printed lines.
"""
import logging
import sys

LOGGER_NAME = 'ecco_insitu'


class _StdoutHandler(logging.StreamHandler):
    """Writes to whatever ``sys.stdout`` is at the time (it changes in tests and when captured)."""

    def emit(self, record):
        self.stream = sys.stdout
        super().emit(record)


def get_logger(name=None):
    """The package logger (or a child, e.g. ``get_logger('pipeline')``), set up on first use:
    level INFO, plain ``%(message)s`` lines on stdout. Nothing is added if the application has
    already configured a handler for it."""
    root = logging.getLogger(LOGGER_NAME)
    if not root.handlers:
        h = _StdoutHandler()
        h.setFormatter(logging.Formatter('%(message)s'))
        root.addHandler(h)
        root.setLevel(logging.INFO)
        root.propagate = False
    return root.getChild(name) if name else root
