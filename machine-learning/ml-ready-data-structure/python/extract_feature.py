'''
extract_feature: Python function for extracting features from signals to create a structure containing all extracted features and relevant info regrading Input_Data and all functions

     Input arguments:
       - input_data_list: List of files in the input data folder and subfolders
       - params: Structure containing parameters for feature extraction
       - data_output_path: Path where the output data will be saved
       - desired_output_file: Name of the desired output file

Date: July 4, 2024
Location: Emory University, Georgia, USA
By: Seyedeh Somayyeh Mousavi
Email: bmemousavi@gmail.com
'''
##########################################################################################################################
##########################################################################################################################
### Part 1: Libraries
import os
import sys
import numpy as np
import scipy.io
import json
import pandas as pd
from scipy.signal import find_peaks
from numpy import log
from scipy.fft import fft, fftfreq
from scipy.signal import welch, find_peaks, hann
from scipy.signal import butter, filtfilt
from numpy.fft import fft, fftshift, fftfreq , ifft
from scipy.signal import find_peaks
from scipy.stats import entropy
from scipy import signal
from scipy.signal import iirnotch, lfilter
from scipy.stats import zscore
import warnings
from collections import defaultdict
from scipy.io import loadmat
from function_callers import *
import shutil
###########################################################################################################################
########################################################################################################################### 
def extract_feature(input_data_list, Input_MD5_list, params, data_output_folder, desired_output_file):

    print("Extract feature for:")
    # Initialize the output variable as a dictionary (Don't change)
    analysis = {}

    # Define the list of methods/functions to apply (Don't change)
    list_functions = params['function_dict']

    # Loop through each file and process the data (Don't change)
    for i, file in enumerate(input_data_list):

        # Get file name and folder (Don't change)
        file_name = os.path.basename(file)
        print(file_name)
        file_path = os.path.dirname(file)  
        # print(file_path)

        # Construct the base filename (without extension)
        basename, _ = os.path.splitext(file_name)

        # Concatenate strings to form the checksum file name
        MD5_Checksum_file_name = f'MD5_{basename}.txt'

        # Search for MD5_Checksum_file_name in the Input_MD5_list['name'] to find the checksum_file_index
        checksum_file_index = Input_MD5_list['name'].tolist().index(MD5_Checksum_file_name) if MD5_Checksum_file_name in Input_MD5_list['name'].tolist() else None

        if checksum_file_index is not None:
            # Get the folder path corresponding to the found index
            MD5_Checksum_file_path = Input_MD5_list.loc[checksum_file_index, 'folder']
      
        # Construct full path to the MD5 checksum file
        MD5_Checksum_file_path = os.path.join(MD5_Checksum_file_path, MD5_Checksum_file_name)
                
        # Check if the MD5 file exists
        if os.path.isfile(MD5_Checksum_file_path):
            # Open the MD5 file and read the checksum
            with open(MD5_Checksum_file_path, "r") as f:
                md5_checksum_file = f.read().strip()      

        # Initialize feature dictionary for this file (Don't change)
        feature = {}

        # Add file info (Don't change)
        feature['name'] = file_name
        feature['rel_path'] = file_path
        feature['md5chsum'] = md5_checksum_file
        feature['sampling_freq'] = params['Fs']
        feature['mains_freq'] = params['Powerline_frequency']
        feature['num_ch'] = params['Number_Channels']
        ###########################################################################################################################  
        # Necessary Function (Don't change)
        # Load the data from the file
        Num_fun = 1
        method = list_functions['Fun_0']
        feature['methods'] = {}
        feature['methods'][f'method_{Num_fun}'] = {
            'name': method['name'],
            'rel_path': method['rel_path'],
            'codebase_md5chsum': method['codebase_md5chsum'],
            'codebase_git_repo': method['codebase_git_repo'],
            'codebase_git_commit_id': method['codebase_git_commit_id']
        }

        # Call the function to process the data
        result_0, data = functioncaller_0(method['name'], file_name, file_path, params)

        # Store the input info of the function
        feature['methods'][f'method_{Num_fun}']['inputs'] = {
            'name': 'whole_data',
            'type': 'double',
            'channels': 0
        }
        '''
        # Store the params of the function
        feature['methods'][f'method_{Num_fun}']['params'] = {
             'name': [],
             'type': [],
             'val': []
         }
                
        # Store the result of the function
        feature['methods'][f'method_{Num_fun}']['outputs'] = {
             'name': [],
             'type': [],
             'val': []
         }
        '''
        feature['methods'][f'method_{Num_fun}']['warning'] = result_0['warning']
        feature['methods'][f'method_{Num_fun}']['success'] = result_0['success']
        ###########################################################################################################################  
        # Part 1: Function_Num1 
        Num_fun += 1
        method = list_functions['Fun_1']
        feature['methods'][f'method_{Num_fun}'] = {
            'name': method['name'],
            'rel_path': method['rel_path'],
            'codebase_md5chsum': method['codebase_md5chsum'],
            'codebase_git_repo': method['codebase_git_repo'],
            'codebase_git_commit_id': method['codebase_git_commit_id']
        }

        channel_name = [None] * params['Number_Channels']
        warning = [None] * params['Number_Channels']
        success = [0] * params['Number_Channels']
        lf_power = [None] * params['Number_Channels']
        hf_power = [None] * params['Number_Channels']
        vlf_power = [None] * params['Number_Channels']
        lf_hf_ratio = [None] * params['Number_Channels']

        # Process each channel individually
        for k in range(params['Number_Channels']):
            
            # Call the second function 
            result_1 = functioncaller_1(method['name'], data, k, file_name, params)
            channel_name[k] = k + 1
            success[k] = result_1['success']
            warning[k] = result_1['warning']
            lf_power[k] = result_1['lf_power']
            hf_power[k] = result_1['hf_power']
            vlf_power[k] = result_1['vlf_power']
            lf_hf_ratio[k] = result_1['lf_hf_ratio']
        
        # Store the input info of the function
        feature['methods'][f'method_{Num_fun}']['inputs'] = {
            'name': 'channel_data',
            'type': 'double',
            'channels': channel_name
        }

        # Store the params of the function
        feature['methods'][f'method_{Num_fun}']['params'] = {
            'name': 'fs',
            'type': 'float',
            'val': params['Fs']
        }
        
        # Store the result of the function
        feature['methods'][f'method_{Num_fun}']['outputs'] = {
            'output_1': {
                'name': 'lf_power',
                'type': 'double',
                'val': lf_power
            },
            'output_2': {
                'name': 'hf_power',
                'type': 'double',
                'val': hf_power
            },
            'output_3': {
                'name': 'vlf_power',
                'type': 'double',
                'val': vlf_power
            },
            'output_4': {
                'name': 'lf_hf_ratio',
                'type': 'double',
                'val': lf_hf_ratio
            }        
            
        }
        feature['methods'][f'method_{Num_fun}']['warning'] = warning
        feature['methods'][f'method_{Num_fun}']['success'] = success 
        ###########################################################################################################################  
        # Part 2: Function_Num2
        Num_fun+=1
        method = list_functions['Fun_2']
        feature['methods'][f'method_{Num_fun}'] = {
                'name': method['name'],
                'rel_path': method['rel_path'],
                'codebase_md5chsum': method['codebase_md5chsum'],
                'codebase_git_repo': method['codebase_git_repo'],
                'codebase_git_commit_id': method['codebase_git_commit_id']
                    }

        # Call the function to process the data
        result_2 = functioncaller_2(method['name'], data, file_name, params)

        # Store the input info of the function
        feature['methods'][f'method_{Num_fun}']['inputs'] = {
            'name': 'whole_data',
            'type': 'double',
            'channels': 0
        }

        # Store the params of the function
        feature['methods'][f'method_{Num_fun}']['params'] = {
             'name': [],
             'type': [],
             'val': []
         }

        # Store the result of the function
        feature['methods'][f'method_{Num_fun}']['outputs'] = {
            'output_1': {
                'name': 'minvalues',
                'type': 'double',
                'val': result_2['minvalues']
                },
            'output_2': {
                'name': 'maxvalues',
                'type': 'double',
                'val': result_2['maxvalues']
                    }
            }        
        # Store the warning and success status
        feature['methods'][f'method_{Num_fun}']['warning'] = result_2['warning']    
        feature['methods'][f'method_{Num_fun}']['success'] = result_2['success'] 
        ############################################################################################################################## 
        # Part 3: Function_Num3
        
        ###############################################################################################################################
        # Store the feature data in the analysis struct (Don't change)
        analysis[f'record_{i + 1}'] = feature

###############################################################################################################################
    # Last part: Save (Don't change)
    # Convert the extracted analysis to JSON format
    json_data = json.dumps(analysis, indent=4)
    
    # Create full path for the output JSON file 
    output_file_path = os.path.join(data_output_folder, desired_output_file)
    
    # Write JSON data to file
    with open(output_file_path, 'w') as f:
            f.write(json_data)

###############################################################################################################################
    # If MD5_Checksum_records exists in the output folder, remove it (Don't change)
    full_path = os.path.join(data_output_folder, "MD5_Checksum_records_python")

    if os.path.isdir(full_path):
       shutil.rmtree(full_path)
