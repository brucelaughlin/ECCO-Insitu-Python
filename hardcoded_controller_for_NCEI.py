from pathlib import Path
import NCEI

output_base_directory = "/Users/brucel/ecco/yip/ECCO-Insitu-Python/processed_profile_files"
#output_base_directory = "/Users/brucel/ecco/yip/sample_data/test_output_all_files"

base_directory = "/Users/brucel/ecco/yip/scripps_data"
directories = [d for d in Path(base_directory).rglob('*') if d.is_dir()]
for input_directory in directories:
    if len([f for f in input_directory.glob('*.nc') if f.is_file()]) > 0:
        output_directory = Path(output_base_directory) / input_directory.name
        output_directory.mkdir(parents=True, exist_ok=True)
        NCEI.NCEI_pipeline(str(output_directory), str(input_directory))
