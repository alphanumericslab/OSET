'''
functioncaller: Includes all Python functions needed to extract features from signals.

Date: July 4, 2024
Location: Emory University, Georgia, USA
By: Seyedeh Somayyeh Mousavi
Email: bmemousavi@gmail.com
'''
###########################################################################################################################
# #########################################################################################################################
# ## Part 1: Libraries
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
###########################################################################################################################
# ## Part 2: Load External Functions
'''
module_path = '../../external/antropy/'
if module_path not in sys.path:
    sys.path.append(module_path)
'''    
###########################################################################################################################
# #########################################################################################################################
def functioncaller_0(method_name, file_name, file_path, params):
    """
    functioncaller_0 - Calls a specified function to load data from a file.

    Inputs:
        method_name: Name of the method/function to call for loading data
        file_name: Name of the file to load data from
        file_path: Path where the file is located
        params: Dictionary containing various parameters including sampling
                frequency, number of channels, and all input parameters of the feature extraction functions

    Outputs:
        feature: Dictionary indicating success of the function call (success = True/False)
        data: Loaded data from the file after processing (if successful)
    
    """
    
    feature = {'success': 0, 'warning': []}
    data = []
    
    try:
        # call the function (Change)
        data = loadmat(os.path.join(file_path, file_name)) # 'pd.read_csv'
        # Modify this part
        data = data['val']  
        # Don't change
        feature['success'] = 1

    except Exception as e:
        # Handle errors during data loading
        feature['warning'] = f"Failed to perform Func: {method_name} for file: {file_name}. Error: {str(e)}"

    # Process data if loading was successful
    if feature['success']:
        # Check the shape of the data and reshape if necessary
        if data.shape[0] != params['Number_Channels']:
            data = data.T  # Transpose data if number of rows doesn't match Number_Channels

    return feature, data

# #########################################################################################################################
def functioncaller_1(method_name, data, k, file_name, Params):
    """
    functioncaller_1 - Calls a specified function to detect peaks in data.

    Inputs:
        method_name: Name of the method/function to call for peak detection
        data: Data on which peak detection is performed
        k: Index of channel
        file_name: Name of the file (used for error reporting)
        Params: Dictionary containing various parameters including sampling
                frequency, number of channels, and all input parameters of the feature extraction functions

    Output:
        channel_feature: Dictionary containing peak detection results and success status
    """
    channel_feature = {'success': 0, 'warning': [], 'lf_power': [], 'hf_power': [], 'vlf_power': [], 'lf_hf_ratio': []}

    try:
        data_channel = data[k, :]
        # call the function
        method = eval(method_name)
        # Call the function to detect peaks
        lf_power, hf_power, vlf_power, lf_hf_ratio = method(data_channel, Params['Fs'])

        # Store results and indicate success
        channel_feature['success'] = 1
        channel_feature['lf_power'] = lf_power
        channel_feature['hf_power'] = hf_power
        channel_feature['vlf_power'] = vlf_power
        channel_feature['lf_hf_ratio'] = lf_hf_ratio

    except Exception as e:
        # Handle errors during peak detection
        channel_feature['warning'] = f"Failed to perform Func: {method_name} for file: {file_name}. Error: {str(e)}"

    return channel_feature
    
# #########################################################################################################################
def functioncaller_2(method_name, data, file_name, Params):
    """
    functioncaller_2 - Calls a specified function to compute extremum values (min and max) from data.

    Inputs:
        method_name: Name of the method/function to call for computing extremum values
        data: 2D numpy array containing data for which extremum values are computed
        file_name: Name of the file (used for error reporting)
        Params: Dictionary containing various parameters including sampling
                frequency, number of channels, and all input parameters of the feature extraction functions

    Outputs:
        feature: Dictionary indicating success of the function call and containing extremum values
    """
    feature = {'success': 0, 'warning': [], 'minvalues': [], 'maxvalues': []}
    
    try:
        # Call the function to compute extremum values
        method = eval(method_name)
        minvalues, maxvalues = method(data)
        
        # Store results and indicate success
        feature['minvalues'] = minvalues
        feature['maxvalues'] = maxvalues
        feature['success'] = 1

    except Exception as e:
        # Handle errors during function execution
        feature['warning'] = f"Failed to perform Func: {method_name} for file: {file_name}. Error: {str(e)}"

    return feature
    
# #########################################################################################################################
def frequency_domain_features(x, Fs):
    """
    frequency_domain_features: Function to calculate frequency domain features from a signal
    
    Inputs:
        x - A 1D numpy array representing the input signal
        Fs - Sampling frequency of the input signal
    
    Outputs:
        lf_power - Power in the low frequency range (0.04-0.15 Hz)
        hf_power - Power in the high frequency range (0.15-0.4 Hz)
        vlf_power - Power in the very low frequency range (0.0033-0.04 Hz)
        lf_hf_ratio - Ratio of low frequency power to high frequency power
    """
    # Shifted FFT
    X_fft_shifted = fftshift(fft(x))
    N = len(X_fft_shifted)
    n = np.arange(N)
    T = 1 / Fs 
    freq_shifted = fftshift(fftfreq(N, T))

    # Power Spectral Density (PSD) using FFT
    psd_fft = np.abs(X_fft_shifted)**2

    # Define frequency ranges for LF, HF, and VLF
    lf_range = (0.04, 0.15)
    hf_range = (0.15, 0.4)
    vlf_range = (0.0033, 0.04)

    # Extract LF, HF, and VLF components using PSD (FFT)
    lf_power = np.sum(psd_fft[(freq_shifted >= lf_range[0]) & (freq_shifted <= lf_range[1])])
    hf_power = np.sum(psd_fft[(freq_shifted >= hf_range[0]) & (freq_shifted <= hf_range[1])])
    vlf_power = np.sum(psd_fft[(freq_shifted >= vlf_range[0]) & (freq_shifted <= vlf_range[1])])

    # Round the values to 3 decimal places
    lf_power = round(lf_power, 3)
    hf_power = round(hf_power, 3)
    vlf_power = round(vlf_power, 3)

    # Calculate LF/HF ratio using PSD (FFT)
    lf_hf_ratio = lf_power / hf_power
    lf_hf_ratio = round(lf_hf_ratio, 3)
    
    return lf_power, hf_power, vlf_power, lf_hf_ratio

# #########################################################################################################################
def Extremum_values(data):
    """
    Extremum_values: Function to find the minimum and maximum values in each column of the input data matrix
    
    Inputs:
        data - A 2D numpy array where we want to find the minimum and maximum values in each column
        
    Outputs:
        minvalues - A list containing the minimum value in each column of the input data
        maxvalues - A list containing the maximum value in each column of the input data
    """
    # Find the minimum values along the first dimension (rows)
    minvalues = np.min(data, axis=0).tolist()
    
    # Find the maximum values along the first dimension (rows)
    maxvalues = np.max(data, axis=0).tolist()
    
    return minvalues, maxvalues

