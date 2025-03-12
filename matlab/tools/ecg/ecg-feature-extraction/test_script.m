% A test script for extracting features from ECG signals in the format of WFDB files
% Date: Jan 10, 2025
% Author: Seyedeh Somayyeh Mousavi
% Location: Emory University, Georgia, USA
% Email: bmemousavi@gmail.com
%===================================================================================================
clc
clear
close all
%% ==================================================================================================
% path
input_path = '../../../sample-data/';
output_path = pwd;

% Name of the output .csv file to save the results
output_file_name = 'ECG_features.csv';

% Define the list of leads
lead_list = {'I', 'II', 'III', 'aVR', 'aVL', 'aVF', 'V1', 'V2', 'V3', 'V4', 'V5', 'V6'};

% Define the number of SVD features
n_eigenvalues = 15;

% Define the notch filter frequency for noise removal
f_notch = 60;

% Run feature extraction for a folder containing WFDB files
extract_ecg_features_path_wfdb(input_path, output_path, output_file_name, lead_list, n_eigenvalues, f_notch)