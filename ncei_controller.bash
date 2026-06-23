

base_input_directory="/Users/brucel/ecco/yip/scripps_data"
input_directory_leaves=($(ls "$base_input_directory"))



output_directory="/Users/brucel/ecco/yip/ECCO-Insitu-Python/processed_profile_files"
mkdir -p "$output_directory"


log_dir="/Users/brucel/ecco/yip/ECCO-Insitu-Python/z_logs"
mkdir -p "$log_dir"


count=1
for input_directory_leaf in ${input_directory_leaves[@]}; do
    #echo "$input_directory"
    #echo $count
    #((count++))

    log_stdout="$log_dir/${input_directory_leaf}_stdout.log"
    log_stderr="$log_dir/${input_directory_leaf}_stderr.log"
    
    input_directory="${base_input_directory}/${input_directory_leaf}"
    
    #echo ""
    #echo $input_directory
    #echo $log_stdout
    #echo $log_stderr


    python -u NCEI.py --input_dir "$input_directory" --dest_dir "$output_directory" 1>"$log_stdout" 2>"$log_stderr" &

done
