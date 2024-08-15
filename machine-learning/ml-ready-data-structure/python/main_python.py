# %%
"""
# General template of Python script for extracting features using various functions from signals 
"""

# %%
"""
### Date: June 19, 2024
### Location: Emory University, Georgia, USA
### By: Seyedeh Somayyeh Mousavi
### Email: bmemousavi@gmail.com
"""

# %%
import os
import json
import hashlib
import scipy.io
from scipy.io import loadmat
from extract_feature import *

# %%
"""
# Part 1: Setup and Initialization (Modify)
"""

# %%
# Part 2.1: Path containing the data (Change)
Data_Input_Path = '../sample_data/'

# Part 2.2: Path where the results will be saved (Change)
Data_Output_Path = '../results/'

# Part 2.3: Name of the output file (Change)
Desired_Output_File = 'Example_ECG_ML_Ready_Python.json'

# Part 2.4: List of files in the folder and subfolders (Change)
data_endswith = '.mat'

# Don't change
Input_Data_list = []
for root, dirs, files in os.walk(Data_Input_Path):
    for file in files:
        if file.endswith(data_endswith):
            Input_Data_list.append(os.path.join(root, file))

# %%
"""
# Part 2: Parameter Configuration (Modify)
"""

# %%
Params = {}

# Part 1.1: Sampling frequency (Hz) (Change)
Params['Fs'] = 500.0

# Part 1.2: Number of channels in the dataset (Change)
Params['Number_Channels'] = 12

# Part 1.3: Powerline frequency (Hz) (Change)
Params['Powerline_frequency'] = 60

# Part 1.4: Define the function details that you will use in the next steps
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
        'name': 'frequency_domain_features',
        'rel_path': '../',
        'codebase_git_repo': '',
        'codebase_git_commit_id': ''
    },
    # Fun_2 (Change)
    'Fun_2': {
        'name': 'Extremum_values',
        'rel_path': '../',
        'codebase_git_repo': '',
        'codebase_git_commit_id': ''
    }
}

# %%
"""
# Part 3: Validations (Don't Change)
"""

# %%
# Check if there are no files in the input data list
if not Input_Data_list:
    raise ValueError('No signal files found in the Data_Input_Folder. Please check the folder and try again.')

# Check if Fs is a numerical value
if not isinstance(Params['Fs'], (int, float)):
    raise ValueError('Fs must be a numerical value.')

# Check if Number_channels is a numerical value 
if not isinstance(Params['Number_Channels'], (int, float)):
    raise ValueError('Number of channels must be a numerical value.')

# Check if the Output_Path exists
if not os.path.isdir(Data_Output_Path):
    raise NotADirectoryError(f'Output_Path "{Data_Output_Path}" does not exist. Please check the folder path and try again.')

# %%
"""
# Part 4.1: Generate_md5 checksums for functions (Don't Change)
"""

# %%
print("Generate_md5 checksums")

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
        print(f"Warning: This path is not a directory: {path}.")

    func_details['codebase_md5chsum'] = checksums

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

    results_folder_path = results_folder_path + "MD5_Checksum_records_python"
    
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

# List of MD5 checksum
Input_MD5_list = []

Input_MD5 = Data_Output_Path + "MD5_Checksum_records_python"

# Walk through the directory tree starting from Data_Input_Path
for root, dirs, files in os.walk(Input_MD5):
    for file in files:
        # Check if the file starts with 'MD5_' and ends with '.txt'
        if file.startswith('MD5_') and file.endswith('.txt'):
            
            # Get the full file path
            file_path = os.path.join(root, file)
            # Get the file stats
            file_stats = os.stat(file_path)
            
            # Append a dictionary to the list
            Input_MD5_list.append({
                'name': file,
                'folder': root,
            })
Input_MD5_list = pd.DataFrame(Input_MD5_list)           

# %%
"""
# Part 5: Feature Extraction (Don't Change)
"""

# %%
"""
## Call the feature extraction function with the specified parameters
#### - Input_Data_list: List of files to process
#### - Params: Struct containing various parameters including sampling frequency, number of channels, and all input variables of the feature extraction functions
#### - Data_Output_Path: Folder to save the results
#### - Desired_Output_File: Name of the output file
"""

# %%
extract_feature(Input_Data_list, Input_MD5_list, Params, Data_Output_Path, Desired_Output_File)