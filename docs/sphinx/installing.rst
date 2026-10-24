Installing
==========

The package is installed from the repository checkout:

.. code-block:: bash

   git clone https://github.com/brucelaughlin/ECCO-Insitu-Python
   cd ECCO-Insitu-Python
   pip install -e .                       # -e: linked to the checkout, so code edits take effect at once
   pip install -e ".[diagnostics,test,docs]"   # also the map/report extras, pytest, and Sphinx

This gives the ``ecco-insitu`` command (``ecco-insitu --help`` lists everything). The code is in
``src/ecco_insitu/``.

**Dependencies** (from ``pyproject.toml``): numpy, xarray, netCDF4, scipy, pandas, pyyaml, pyproj,
geopy, pymatreader, matplotlib. Extras: ``diagnostics`` (cartopy, pyarrow), ``test`` (pytest),
``docs`` (sphinx, sphinx-rtd-theme, myst-parser).

**Data the chain needs** (not in the repository; a support package is on the project Google Drive):
the LLC90 grid, the geodesic bin files, the CTD-based T/S uncertainty fields, the Arctic
high-latitude uncertainty files, and a climatology (which can be built, see :doc:`quickstart`).
Their locations go in ``config/local.yaml`` (:doc:`settings`).

**Tests** (about 75 s; two of them run the whole chain on small test files):

.. code-block:: bash

   python -m pytest
