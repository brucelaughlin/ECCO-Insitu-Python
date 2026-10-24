"""
Step 08: optionally trim trailing empty depth levels.

Off by default (``chain.remove_extraneous_depth_levels: false``), so every output file stays on
the full input depth grid. When on, the depth levels below the deepest level with valid data, in
any variable of the file, are removed.
"""
# Probably my favorite module ever

import ecco_insitu.tools as tools
def main(MITprof_ds, profile_var_key_set, remove_extraneous_depth_levels=True):
    """Trim trailing depth levels with no valid data in the file.

    Parameters
    ----------
    MITprof_ds : xarray.Dataset
        The profiles.
    profile_var_key_set : set of str
        Variables whose data count (``prof_T``, ``prof_S``).
    remove_extraneous_depth_levels : bool
        Trim (``chain.remove_extraneous_depth_levels``; false by default).

    Returns
    -------
    xarray.Dataset
        Unchanged when not trimming; otherwise without the levels below the deepest one with
        data.

    Notes
    -----
    When not trimming, every output file stays on the full input depth grid: empty levels are
    written as fill, and their weights are already zero from step 06. The step stays in the chain
    either way, so the step numbering (step 09 = decimation, file names, log messages) doesn't
    change.
    """
    if not remove_extraneous_depth_levels:
        print("  step 08 disabled (remove_extraneous_depth_levels = False): keeping all depth levels")
        return MITprof_ds
    #print("step08: update_remove_extraneous_depth_levels")
    MITprof_ds = tools.update_remove_extraneous_depth_levels(MITprof_ds, profile_var_key_set) 

    return MITprof_ds
