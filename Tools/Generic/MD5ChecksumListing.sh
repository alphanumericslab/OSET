# This script recursively searches a root folder for files with a given extension and logs the MD5 sum of the files and their full path in a CSV file
#!/bin/bash

root_folder="/labs/samenilab/data/physionet.org/files/challenge-2017/1.0.0/training/A01/"
output_file="./log_file.csv"
# Example:
#   root_folder="/Users/rsameni/Documents/DataFiles/physionet.org/files/challenge-2017/1.0.0/training/A00"
#   output_file="./log_file.csv"

file_extension=".mat" # file extensions of interest

os=$(uname -s)
if [[ "$os" == "Linux" ]]; then # FOR LINUX:
    # Recursively search for all files with the specified extension in root_folder and calculate their MD5 sums
    find "$root_folder" -type f -name "*$file_extension" -exec md5sum {} \; | awk '{gsub(/ /,",",$0); sub(/,/,"",$0); print}' > "$output_file"

elif [[ "$os" == "Darwin" ]]; then # FOR MAC OS
    # Recursively search for all files with the specified extension in root_folder and calculate their MD5 sums
    find "$root_folder" -type f -name "*$file_extension" -exec md5 -r {} \; > "$output_file"

    # Replace spaces with commas in the output file to create a CSV file
    sed -i '' 's/ /,/g' "$output_file"
else
    echo "This code only runs on Linux or macOS."
fi
