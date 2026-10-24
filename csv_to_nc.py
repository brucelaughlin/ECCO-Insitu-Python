"""Compatibility wrapper: this script moved into the ecco_insitu package.

Same as:  ecco-insitu csv-to-nc [options]
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ecco_insitu.cli import main  # noqa: E402

if __name__ == '__main__':
    sys.exit(main(['csv-to-nc', *sys.argv[1:]]))
