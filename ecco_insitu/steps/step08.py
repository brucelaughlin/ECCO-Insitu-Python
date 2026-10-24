# Probably my favorite module ever

import ecco_insitu.tools as tools
def main(MITprof_ds, profile_var_key_set, remove_extraneous_depth_levels=True):
    """Trim trailing depth levels with no valid data in the file.

    Controlled by remove_extraneous_depth_levels in NCEI.py.  When False, the
    step passes the data through untouched so every output file stays on the
    full input depth grid (empty levels are written as fill, and their weights
    are already zero from step 06).  The step stays in the chain either way so
    step numbering (step 09 = decimation, file names, log messages) is unchanged.
    """
    if not remove_extraneous_depth_levels:
        print("  step 08 disabled (remove_extraneous_depth_levels = False): keeping all depth levels")
        return MITprof_ds
    #print("step08: update_remove_extraneous_depth_levels")
    MITprof_ds = tools.update_remove_extraneous_depth_levels(MITprof_ds, profile_var_key_set) 

    return MITprof_ds
