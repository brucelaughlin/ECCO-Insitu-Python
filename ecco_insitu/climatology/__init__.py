"""The climatology the chain compares profiles with (step 03): check, build, pre-bake.

You don't need a climatology already. Either point ``paths.climatology_file`` at your own file
(the contract is in ``docs/input_contract.md``, section 2), or build the standard WOA23 one with
``ecco-insitu climatology --build``. See :mod:`ecco_insitu.climatology.setup`.
"""
