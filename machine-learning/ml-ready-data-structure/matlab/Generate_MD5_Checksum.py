# %%
"""
# Python script for Generating_md5 checksums
"""

# %%
"""
### Date: Junly 29, 2024
### Location: Emory University, Georgia, USA
### By: Seyedeh Somayyeh Mousavi
### Email: bmemousavi@gmail.com
"""

# %%
import os
import json
import hashlib
import scipy.io
import pandas as pd
from scipy.io import loadmat

# %%
"""
# Part 1: Setup and Initialization (Modify)
"""

# %%
# Part 1.1: Path containing the data (Change)
Data_Input_Path = '../sample_data/'

# Part 1.2: Path where the results will be saved (Change)
Data_Output_Path = '../results/'

# Part 1.3: Create the folder for MD5 Checksum records (Don't change)
MD5_Checksum_Path = os.path.join(Data_Output_Path, 'MD5_Checksum_records_matlab') 
# Create the directory if it doesn't exist (Don't change)
os.makedirs(MD5_Checksum_Path, exist_ok=True)

# Part 1.4: List of files in the folder and subfolders (Change)
data_endswith = '.mat' 
# (Don't change)
Input_Data_list = []
for root, dirs, files in os.walk(Data_Input_Path):
    for file in files:
        if file.endswith(data_endswith):
            Input_Data_list.append(os.path.join(root, file))

# %%
"""
# Part 2: Validations (Don't Change)
"""

# %%
# Check if there are no files in the input data list
if not Input_Data_list:
    raise ValueError('No signal files found in the Data_Input_Folder. Please check the folder and try again.')

# Check if the Output_Path exists
if not os.path.isdir(Data_Output_Path):
    raise NotADirectoryError(f'Output_Path "{Data_Output_Path}" does not exist. Please check the folder path and try again.')

# %%
"""
# Part 3: Parameter Configuration (Modify)
"""

# %%
Params = {}

# Define the function details that you will use in the next steps
Params['function_dict'] = {
    
    # Fun_0: (Don't Change)
    'Fun_0': {
        'name': 'load', 
        'rel_path': '',
        'codebase_git_repo': '',
        'codebase_git_commit_id': ''
    },
    # Fun_1 (Change)
    'Fun_1': {
        'name': 'peak_det_modified_pan_tompkins',
        'rel_path': '../',
        'codebase_git_repo': '',
        'codebase_git_commit_id': ''
    },
    # Fun_2 (Change)
    'Fun_2': {
        'name': 'max_values',
        'rel_path': '../',
        'codebase_git_repo': '',
        'codebase_git_commit_id': ''
    }
}

# %%
"""
# Part 4.1: Generate_md5 checksums for functions (Don't Change)
"""

# %%
# Function to calculate the MD5 checksum of a file
def calculate_md5_for_file(file_path):
    md5_hash = hashlib.md5()
    with open(file_path, "rb") as f:
        for byte_block in iter(lambda: f.read(4096), b""):
            md5_hash.update(byte_block)
    return md5_hash.hexdigest()

# Function to calculate a combined MD5 checksum for all files in a directory
def calculate_md5_for_folder(folder_path):
    md5_hash = hashlib.md5()
    for root, dirs, files in os.walk(folder_path):
        # Sorting to ensure consistent order
        for file_name in sorted(files):  
            file_path = os.path.join(root, file_name)
            with open(file_path, "rb") as f:
                for byte_block in iter(lambda: f.read(4096), b""):
                    md5_hash.update(byte_block)
    return md5_hash.hexdigest()

# Add MD5 checksum for each function based on rel_path
for key, func_details in Params['function_dict'].items():
    checksums = ''
    path = func_details['rel_path']

    if os.path.isdir(path):
        # If the path is a directory, calculate MD5 for the whole directory
        checksums = calculate_md5_for_folder(path)
    elif os.path.isfile(path):
        # If the path is a file, calculate MD5 for the file
        checksums = calculate_md5_for_file(path)
    else:
        # Handle the case where the path is neither a file nor a directory
        print(f"Warning: Path is not a directory, {path}.")
        
    func_details['codebase_md5chsum'] = checksums

# %%
"""
# Save the Params (Don't Change)
"""

# %%
# Convert the Params['function_dict'] dictionary to a list of dictionaries
data = []
for key, func_details in Params['function_dict'].items():
    row = {
        'function_key': key,
        'name': func_details['name'],
        'rel_path': func_details['rel_path'],
        'codebase_git_repo': func_details['codebase_git_repo'],
        'codebase_git_commit_id': func_details['codebase_git_commit_id'],
        'codebase_md5chsum': func_details['codebase_md5chsum']
    }
    data.append(row)

# Create a DataFrame from the list of dictionaries
df = pd.DataFrame(data)

# Save the DataFrame as a CSV file
csv_path = os.path.join(MD5_Checksum_Path, "Params_function_dict.csv")
df.to_csv(csv_path, index=False)
print(f"Params_function_dict has been saved as a CSV file to {csv_path}")

# %%
"""
# Part 4.2: Generate_md5 checksums for records (Don't Change)
"""

# %%
def generate_md5_for_record(file_path, output_path):
    
    # Create an MD5 hash object
    hash_md5 = hashlib.md5()
    
    # Open the file in binary read mode
    with open(file_path, "rb") as f:
        
        # Read the file in chunks of 4096 bytes
        for chunk in iter(lambda: f.read(4096), b""):
            # Update the MD5 hash object with the chunk
            hash_md5.update(chunk)
    
    # Get the hexadecimal representation of the MD5 hash
    md5_hash = hash_md5.hexdigest()
    
    # Open the output file in write mode
    with open(output_path, "w") as out_file:
        # Write the MD5 hash to the output file
        out_file.write(md5_hash)
    
    # Return the MD5 hash
    return md5_hash

def generate_md5_for_all_records_in_folder(input_folder_path, results_folder_path):

    results_folder_path = results_folder_path + "MD5_Checksum_records_matlab"
    
    # Walk through the directory tree starting from input_folder_path
    for root, dirs, files in os.walk(input_folder_path):
        
        # Compute the relative path from the input folder path
        relative_path = os.path.relpath(root, input_folder_path)
        
        # Create corresponding folder in the results directory
        result_folder = os.path.join(results_folder_path, relative_path)
        os.makedirs(result_folder, exist_ok=True)
        
        # Iterate over all files in the current directory
        for filename in files:
            
            # Check if the file has the desierd extension
            if filename.endswith(data_endswith):
                
                # Construct the full file path
                file_path = os.path.join(root, filename)
                
                # Construct the output file path
                output_path = os.path.join(result_folder, f"MD5_{os.path.splitext(filename)[0]}.txt")
                
                # Print the filename being processed
                print(f"Generating MD5 for {filename}")
                
                # Generate the MD5 hash for the file and save it
                generate_md5_for_record(file_path, output_path)

generate_md5_for_all_records_in_folder(Data_Input_Path, Data_Output_Path)