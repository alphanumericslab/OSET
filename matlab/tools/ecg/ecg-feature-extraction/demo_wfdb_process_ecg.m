% Description: Demo script for processing ECG data from WFDB files
%
% This script demonstrates how to use the process_ecg_wfdb function to process
% different types of ECG recordings in various formats (.mat, .dat) and save
% the extracted features to CSV files.
%
% The script includes three examples:
% 1. Processing a .mat wfdb file from PTB-XL
% 2. Processing a .dat wfdb file with two leads from QTDB
% 3. Processing a .dat wfdb file from PTB 
% Author: Sajjad Karimi
% Date: Apr 10, 2025
% Location: Emory University, Georgia, USA
% Email: sajjadkarimi91@gmail.com

clc
clear
close all

% Example 1: Process a .mat file with standard 12-lead ECG
% This example processes a short segment (10 seconds) of a 12-lead ECG recording
% Features are stacked by channel (flatten_flag = 0)
input_wfdb_address = './sample-ecg/HR00001.mat';
ecg_csv_file_name = './sample-ecg/HR00001_ecg.csv';
lead_names_target = {'I', 'II', 'III', 'aVR', 'aVL', 'aVF', 'V1', 'V2', 'V3', 'V4', 'V5', 'V6'};
window_dur = 10;    % Process in 10-second windows
start_time = 1;     % Start from 1 second
stop_time = 0;      % Process until the end
flatten_flag = 0;   % Stack features by channel

process_ecg_wfdb(input_wfdb_address, ecg_csv_file_name, lead_names_target, window_dur, start_time, stop_time, flatten_flag);

% Example 2: Process a .dat file with custom lead configuration
% This example processes a longer segment (7 minutes) of a recording with
% custom lead configuration. Features are flattened across channels.
input_wfdb_address = './sample-ecg/sel30.dat';
ecg_csv_file_name = './sample-ecg/sel30_ecg.csv';
lead_names_target = {'ECG1', 'ECG2'};  % Custom lead names
window_dur = 60;    % Process in 1-minute windows
start_time = 10;    % Start from 10 seconds
stop_time = 430;    % Process until 430 seconds
flatten_flag = 1;   % Flatten features across channels

process_ecg_wfdb(input_wfdb_address, ecg_csv_file_name, lead_names_target, window_dur, start_time, stop_time, flatten_flag);

% Example 3: Process a .dat file with reduced lead set
% This example processes a short segment of a recording with a reduced
% lead set. Features are stacked by channel.
input_wfdb_address = './sample-ecg/s0010_re.dat';
ecg_csv_file_name = './sample-ecg/s0010_re_ecg.csv';
lead_names_target = {'i', 'ii', 'ii', 'v1', 'v2', 'v3', 'v4', 'v5', 'v6'};
window_dur = 10;    % Process in 10-second windows
start_time = 0;     % Start from the beginning
stop_time = 0;      % Process until the end
flatten_flag = 0;   % Stack features by channel

process_ecg_wfdb(input_wfdb_address, ecg_csv_file_name, lead_names_target, window_dur, start_time, stop_time, flatten_flag);

% NOTES:
% 1. The flatten_flag parameter affects how features are organized in the output:
%    - flatten_flag = 0: Features are stacked by channel [C x N]
%    - flatten_flag = 1: Features are concatenated across channels [1 x (C*N)]
%
% 2. Window duration (window_dur) affects processing time and memory usage:
%
% 3. Lead names must match those in the WFDB header file:
%    - Case sensitive
%
% 4. Time parameters (start_time, stop_time) are in seconds:
%    - 0 means start/end of the recording
%    - Must be within the recording duration
%
% 5. The script assumes the WFDB files are in the ./sample-ecg/ directory:
%    - Modify the paths if your files are in a different location
%    - Ensure both .dat and .hea files are present for WFDB format