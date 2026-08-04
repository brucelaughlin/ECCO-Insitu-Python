import argparse
import glob
import os
import numpy as np
import numpy.ma as ma
import tools 



def main(MITprofs):
    #print("step09: update_remove_extraneous_depth_levels")
    tools.update_remove_extraneous_depth_levels(MITprofs) 

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    
    parser.add_argument("-m", "--MIT_dir", action= "store",
                    help = "File path to NETCDF files containing MITprofs info." , dest= "MIT_dir",
                    type = str, required= True)

    args = parser.parse_args()
    MITprofs_fp = args.MIT_dir
    
    nc_files = glob.glob(os.path.join(MITprofs_fp, '*.nc'))
    if len(nc_files) == 0:
        raise Exception("Invalid NC filepath")
    for file in nc_files:
        MITprofs = MITprof_read(file, 9)

    # Convert all masked arrs to non-masked types
    for keys in MITprofs.keys():
        if ma.isMaskedArray(MITprofs[keys]):
            MITprofs[keys] = MITprofs[keys].filled(np.NaN)
    
    main(MITprofs)
