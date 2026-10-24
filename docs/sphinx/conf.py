"""Sphinx configuration for the ecco_insitu documentation (built from docs/sphinx/)."""
import os
import sys

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
sys.path.insert(0, os.path.join(ROOT, 'src'))

import ecco_insitu  # noqa: E402

project = 'ecco_insitu'
author = 'Sweet Zhang, Ian Fenty, Bruce Laughlin'
copyright = '2023-2026, ' + author
release = ecco_insitu.__version__
version = release

extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.autosummary',
    'sphinx.ext.napoleon',       # NumPy-style docstrings
    'sphinx.ext.viewcode',       # "[source]" links to the highlighted code
    'sphinx.ext.intersphinx',
    'myst_parser',               # the Markdown pages (input contract, diagnostics, ...)
]
autosummary_generate = True
autodoc_default_options = {'members': True, 'undoc-members': False, 'show-inheritance': True}
autodoc_member_order = 'bysource'
# The diagnostics modules import cartopy and pyarrow (the 'diagnostics' extra): install it to build
# the docs (pip install -e ".[diagnostics,docs]"). They are not mocked: a mocked cartopy breaks
# xarray's backend discovery on import.
napoleon_numpy_docstring = True
napoleon_google_docstring = False
# Document a class's attributes from its 'Attributes' docstring section only (dataclass fields would
# otherwise be listed twice).
napoleon_use_ivar = True
myst_heading_anchors = 3
source_suffix = {'.rst': 'restructuredtext', '.md': 'markdown'}
intersphinx_mapping = {
    'python': ('https://docs.python.org/3', None),
    'numpy': ('https://numpy.org/doc/stable/', None),
    'xarray': ('https://docs.xarray.dev/en/stable/', None),
}
html_theme = 'sphinx_rtd_theme'
html_title = f'ecco_insitu {release}'
exclude_patterns = ['_build']
templates_path = ['_templates']
