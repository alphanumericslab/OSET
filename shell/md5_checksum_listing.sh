#!/bin/bash

# Script to recursively search a root folder for files with a given extension and log their MD5 sum and relative path in a CSV file
# Usage: md5_checksum_listing.sh root_folder output_file file_extension
# Example: bash md5_checksum_listing.sh /labs/samenilab/data/physionet.org/files/challenge-2017/1.0.0/training/ ./RECORDS.csv .mat
# 
# Note: Remember to add write permission chmod +x md5_checksum_listing.sh

search_files_md5sum() {
    local root_folder="$1"
    local output_file="$2"
    local file_extension="$3" # File extension of interest

    if [[ -z "$root_folder" || -z "$output_file" || -z "$file_extension" ]]; then
        echo "Usage: md5_checksum_listing.sh root_folder output_file file_extension"
        return 1
    fi

    os=$(uname -s)

    if [[ "$os" == "Linux" ]]; then # For Linux:
        find "$root_folder" -type f -name "*$file_extension" -exec md5sum {} \; | awk -v root="$root_folder" -v ext="$file_extension" '{gsub(root,"",$2); gsub(ext,"",$2); print $1 "," "\"" "/" $2 "\"" }' > "$output_file"
    elif [[ "$os" == "Darwin" ]]; then # For macOS:
        find "$root_folder" -type f -name "*$file_extension" -exec md5 -r {} \; | awk -v root="$root_folder" -v ext="$file_extension" '{gsub(root,"",$2); gsub(ext,"",$2); print $1 "," "\"" $2 ext "\""}' > "$output_file"
    else
        echo "This code only runs on Linux or macOS."
        return 1
    fi
}

# Check if the script is called with the correct number of arguments or help option
if [[ $# -eq 1 && ($1 == "-h" || $1 == "--help") ]]; then
    echo "md5_checksum_listing.sh is a script to recursively search a root folder for files with a given extension and log their MD5 sum and relative path in a CSV file"
    echo "Usage: md5_checksum_listing.sh root_folder output_file file_extension"
    echo "Example: bash md5_checksum_listing.sh /labs/samenilab/data/physionet.org/files/challenge-2017/1.0.0/training/ ./RECORDS.csv .mat"
    echo "Reza Sameni, 2023"
    echo "The Open-Source Electrophysiological Toolbox"
    echo "https://github.com/alphanumericslab/OSET"

    exit 0
elif [[ $# -ne 3 ]]; then
    echo "Invalid arguments. Use -h or --help for usage instructions."
    exit 1
fi

# Call the search_files_md5sum function with the provided arguments
search_files_md5sum "$1" "$2" "$3"
