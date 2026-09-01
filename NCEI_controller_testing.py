from pathlib import Path
import NCEI

output_base_directory = "/Users/brucel/ecco/yip/sample_data/test_output_problematic_files"

input_directory = "/Users/brucel/ecco/yip/ECCO-Insitu-Python/input_files_should_be_easy"
#input_directory = "/Users/brucel/ecco/yip/ECCO-Insitu-Python/input_files_problematic_test"
#input_directory = "/Users/brucel/ecco/yip/ECCO-Insitu-Python/input_files_biggun"

if len([f for f in Path(input_directory).glob('*.nc') if f.is_file()]) > 0:
    output_directory = Path(output_base_directory) 
    output_directory.mkdir(parents=True, exist_ok=True)
    NCEI.NCEI_pipeline(str(output_directory), str(input_directory))
