from pathlib import Path
import NCEI_test_uncertainty

output_base_directory = "/Users/brucel/ecco/yip/sample_data/test_output_problematic_files"

#input_directory = "/Users/brucel/ecco/yip/ECCO-Insitu-Python/input_files_should_be_easy"
#input_directory = "/Users/brucel/ecco/yip/ECCO-Insitu-Python/input_files_problematic_test"
input_directory = "/Users/brucel/ecco/yip/ECCO-Insitu-Python/input_files_biggun"

output_directory = Path(output_base_directory) 
output_directory.mkdir(parents=True, exist_ok=True)
MITprof_ds = NCEI_test_uncertainty.NCEI_pipeline(str(output_directory), str(input_directory))
