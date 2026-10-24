Settings
========

All paths and options are in YAML files, never in the code. They are merged in this order (later
ones win):

1. ``src/ecco_insitu/defaults.yaml``: **every** setting, with its default and a comment explaining it
   (reproduced below);
2. ``config/local.yaml``: this machine's data locations (not committed; copy
   ``config/local.example.yaml``);
3. the file named by the environment variable ``ECCO_INSITU_CONFIG``, if set;
4. each ``--config FILE`` given on the command line;
5. each ``--set key.sub=value``, e.g. ``--set chain.new_floor.prof_S=0.006``.

An unknown setting, or a value of the wrong type, stops the command with a message naming it.
Relative paths are taken from ``data_root``. Each run saves the merged settings it used to
``<run>/documentation/run_config.yaml``.

To see the merged settings: ``ecco-insitu config``; one value: ``ecco-insitu config --get chain.method``.

The loader is :mod:`ecco_insitu.config`.

All settings (``defaults.yaml``)
--------------------------------

.. literalinclude:: ../../src/ecco_insitu/defaults.yaml
   :language: yaml
