
import sys
from pathlib import Path
geodesic_dir = str(Path(__file__).parent.resolve())
sys.path.append(geodesic_dir)
import geodesic_binning_utilities as utils
import in_progress_qc
#import qc

# Turns out it wasn't the number of profiles freezing the program.  Still, we may want something like this....?  For now, 
# leave it in since I'm lazy, just set the threshold high
num_subpolygons_max = 1000
#num_subpolygons_max = 10
#num_subpolygons_max = 100


variables_of_interest = ["T", "S"]
#angular_precision = 0.1
angular_precision = 0.25
#angular_precision = 0.5

geodesic_file_dict = {
    "00642": "/Users/brucel/ecco/yip/sample_data/ecco-insitu/sweet_gdrive/geodesic/00642_bin_locations.csv",
    "02562": "/Users/brucel/ecco/yip/sample_data/ecco-insitu/sweet_gdrive/geodesic/02562_bin_locations.csv",
    "10242": "/Users/brucel/ecco/yip/sample_data/ecco-insitu/sweet_gdrive/geodesic/10242_bin_locations.csv",
}
num_geodesic_bins = 10242

profile_file_list = [
    "/Users/brucel/ecco/yip/ECCO-Insitu-Python/z_test_output/ITP_WO_2004_CTD__ncei_step_10.nc",
    "/Users/brucel/ecco/yip/ECCO-Insitu-Python/z_test_output/WOD_WO_1992_CTD_OSD__ncei_step_10.nc",
    "/Users/brucel/ecco/yip/ECCO-Insitu-Python/z_test_output/WOD_WO_2002_GLD__ncei_step_10.nc",
    "/Users/brucel/ecco/yip/ECCO-Insitu-Python/z_test_output/ARGO_WO_2001_PFL_A__ncei_step_10.nc",
    "/Users/brucel/ecco/yip/ECCO-Insitu-Python/z_test_output/ARGO_WO_2015_PFL_A__ncei_step_10.nc",
    "/Users/brucel/ecco/yip/ECCO-Insitu-Python/z_test_output/ITP2_WO_2008_CTD__ncei_step_10.nc",
]

'''
Good files:
1,4,2

Bad files:
5

Files to test:
3
0
'''

profile_file_index = 2

geodesic_file = geodesic_file_dict[f"{num_geodesic_bins:05}"]
profile_file = profile_file_list[profile_file_index]

geodesic_bin_data = utils.bin_around_geodesic_vertices(geodesic_file, profile_file, variables_of_interest, angular_precision, num_geodesic_bins)


in_progress_qc.pp(geodesic_bin_data, num_geodesic_bins, geodesic_file, profile_file, num_subpolygons_max)
#qc.pp(geodesic_bin_data, num_geodesic_bins)
