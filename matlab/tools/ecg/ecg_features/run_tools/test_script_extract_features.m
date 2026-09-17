% An example script for extracting features from ECG signals
% Authors: 
% Seyedeh Somayyeh Mousavi
% Reza Sameni
% bmemousavi@gmail.com
% Aug 2026
% Emory University, Georgia, USA

% ====================================================================
clc
clear
close all

% ====================================================================
% Input and output paths
input_path = '../../../../../datasets/sample-data/HR00001/';
output_path = './sample/result/HR00001/';
% ====================================================================
% Run the feature extraction code
disp('Run the feature extraction script.')
extract_ecg_features(input_path, output_path)
disp('Running script finished successfully.')
