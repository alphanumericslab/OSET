#!/bin/bash
# This bash script recursively searches a root folder for files with a given extension and logs the MD5 sum of the files and their relative path in a CSV file

root_folder="/labs/samenilab/data/physionet.org/files/challenge-2017/1.0.0/training/"
#root_folder="/Users/rsameni/Documents/DataFiles/physionet.org/files/challenge-2017/1.0.0/training/"
output_file="./RECORDS.csv"

#file_extension=".hea" # file extensions of interest
file_extension=".mat" # file extensions of interest

os=$(uname -s)
# Recursively search for all files with the specified extension in root_folder and calculate their MD5 sums
if [[ "$os" == "Linux" ]]; then # For Linux:
    # Without file extensions
    # find "$root_folder" -type f -name "*$file_extension" -exec md5sum {} \; | awk -v root="$root_folder" -v ext="$file_extension" '{gsub(root,"",$2); gsub(ext,"",$2); print $1 "," "\"" "/" $2 "\"" }' > "$output_file"

    # With file extensions
    find "$root_folder" -type f -name "*$file_extension" -exec md5sum {} \; | awk -v root="$root_folder" -v ext="$file_extension" '{gsub(root,"",$2); gsub(ext,"",$2); print $1 "," "\"" "/" $2 ext "\"" }' > "$output_file"

elif [[ "$os" == "Darwin" ]]; then # For MacOS:
    # Without file extensions
    # find "$root_folder" -type f -name "*$file_extension" -exec md5 -r {} \; | awk -v root="$root_folder" -v ext="$file_extension" '{gsub(root,"",$2); gsub(ext,"",$2); print $1 "," "\"" $2 "\""}' > "$output_file"

    # With file extensions
    find "$root_folder" -type f -name "*$file_extension" -exec md5 -r {} \; | awk -v root="$root_folder" -v ext="$file_extension" '{gsub(root,"",$2); gsub(ext,"",$2); print $1 "," "\"" $2 ext "\""}' > "$output_file"
else
    echo "This code only runs on Linux or macOS."
fi
