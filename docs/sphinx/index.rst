ecco_insitu: in-situ ocean profiles for ECCO
================================================

``ecco_insitu`` takes in-situ temperature and salinity profiles (WOD CTD, XBT, gliders and moorings,
Argo, MEOP seal tags, ITP) through quality control, comparison with a monthly climatology,
uncertainty assignment and spatial scaling, and packs them into the profile files the ECCO ocean
state estimate (MITgcm ``profiles`` package) reads.

.. important::

   **You don't need a climatology already.** The chain compares every profile with a monthly
   potential-temperature / salinity climatology (``paths.climatology_file``). Either point that
   setting at **your own file** (any file meeting :doc:`the climatology contract <input_contract>`),
   or let the package **build the standard one** this project uses (WOA23 1991–2020, 1°, monthly to
   1500 m stitched onto seasonal fields to 5500 m, converted to potential temperature):

   .. code-block:: bash

      ecco-insitu climatology --build      # downloads WOA23 (~1.4 GB) and builds it
      ecco-insitu climatology              # check the setup

   If the climatology is missing or unusable, a run stops before processing and says so.

.. toctree::
   :maxdepth: 2
   :caption: User guide

   quickstart
   installing
   settings
   running
   input_contract
   outputs
   diagnostics
   steps
   model_comparison

.. toctree::
   :maxdepth: 1
   :caption: Reference

   api/index
   automation
   known_issues
   changelog

Indices: :ref:`genindex` · :ref:`modindex`
