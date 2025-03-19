% A test script for extracting all ECG features and saving in .csv format
% Author: Sajjad Karimi
% Location: Emory University, Georgia, USA
% Email: sajjadkarimi91@gmail.com
% Date: Mar 14, 2025

clc
clear
close all
%%
% path
% load sample PTB data from OSET
load("ptb_s0010_re.mat") % load signals, fs and sig_name
fs = double(fs);
output_path = pwd;

% Define the list of leads
for c = 1:size(sig_name,1)
    lead_names{c} = sig_name(c,:);
end

ecg_data = signals';

%% Preprocessing

[C, T] = size(ecg_data); % # of ECG channels
for c = 1:C

    ecg_raw = ecg_data(c,:);
    % Notch filtering ECG
    fc = 60.0; % powerline frequency
    Qfactor = 45; % Q-factor of the notch filter
    Wo = fc/(fs/2);  BW = Wo/Qfactor; % nothc filter parameters
    [b,a] = iirnotch(Wo, BW); % design the notch filter
    ecg_raw = (filtfilt(b, a, ecg_raw))'; % zero-phase non-causal filtering

    % baseline removing
    ecg_raw = ecg_raw - movmean(movmedian(ecg_raw,[round(0.3*fs),round(0.3*fs)]),[round(0.15*fs),round(0.15*fs)]);

    ecg_raw = ecg_raw - lp_filter_zero_phase(ecg_raw, 0.1/fs); % High pass filter
    ecg_raw = lp_filter_zero_phase(ecg_raw, 30/fs); % Low pass filter

    ecg_data(c,:) = ecg_raw;

end

%% Feature extraction
win_len = 10*fs;

num_windows = ceil(T/win_len);

% Name of the output .csv file to save the results
feature_csv_file_name = './s0010_sample_ecg_features.csv';

% Name of the output .csv file to save the results
fiducial_csv_file_name = '/Users/skari25/Documents/projects/toolboxes/OSET/matlab/tools/ecg/ecg-feature-extraction/s0010_sample_ecg_fiducial.csv';

for n = 1: num_windows

    start_index = (n-1)*win_len+1;
    if n<num_windows
        stop_index = n*win_len;
    else
        stop_index = T;
    end

    win_ecg_data = ecg_data(:,start_index:stop_index);
    flatten_flag = false;
    [ecg_features_vector, ecg_feature_info, ecg_fiducial_position, exit_flag] = ecg_feature_extraction(win_ecg_data, fs, lead_names,  [], [], [], [], flatten_flag);


    features_matrix = ecg_features_vector;
    feature_names = ecg_feature_info.names;
    features_units = ecg_feature_info.units;
    feature_description = ecg_feature_info.description;
    window_time_info = repmat( [start_index,stop_index], size(ecg_data,1), 1);

    feature_handle_csv(feature_csv_file_name, features_matrix, feature_names, features_units, feature_description,  fs, window_time_info, lead_names(:))

    ecg_fiducial_handle_csv(fiducial_csv_file_name, win_ecg_data, fs,  ecg_fiducial_position, lead_names, start_index)

end

