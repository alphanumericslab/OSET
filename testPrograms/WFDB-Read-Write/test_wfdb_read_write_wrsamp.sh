
raw_data_file="sample_raw.csv"
channel_gains_file="sample_gain.csv"

# Read CSV file and generate the string with truncated numbers
gain_list_string=$(awk 'BEGIN{FS=","; OFS=" "; ORS="";} {for(i=1;i<=NF;i++){ printf "%.4f", $i; if(i<NF) printf " " } printf "\n" }' "$channel_gains_file")

# Add double quotes around the string
gain_list_string="\"$gain_list_string\""

# echo "$gain_list_string"

# Build the command line and execute it
wrsamp_data_file=$(echo "$raw_data_file" | sed 's/_raw.csv/_wrsamp/')
#echo "$wrsamp_data_file"
cmd="wrsamp -i $raw_data_file -o $wrsamp_data_file -G $gain_list_string -x $gain_list_string -F 256 -O 16"
#echo "$cmd"
eval "$cmd"

# Build the command line and execute it
wrsampNG_data_file=$(echo "$raw_data_file" | sed 's/_raw.csv/_wrsampNG/')
cmd="wrsamp -i $raw_data_file -o $wrsampNG_data_file -G 1 -F 256 -O 16"
#echo "$cmd"
eval "$cmd"
