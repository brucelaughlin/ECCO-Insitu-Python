# Probably my favorite module ever

import tools 

#def main(MITprof_ds, profile_var_key_set=None):
def main(MITprof_ds, profile_var_key_set):
    #print("step09: update_remove_extraneous_depth_levels")
    MITprof_ds = tools.update_remove_extraneous_depth_levels(MITprof_ds, profile_var_key_set) 

    return MITprof_ds
