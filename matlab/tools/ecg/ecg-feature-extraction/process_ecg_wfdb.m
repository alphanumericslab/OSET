function process_ecg_wfdb(input_wfdb_address, ecg_csv_file_name, lead_names_target, window_dur, start_time, stop_time, flatten_flag)

% Description: Process WFDB files to extract ECG features and save them as CSV files
%
% This function reads ECG data from WFDB files, processes it in windows,
% extracts features using ecg_feature_extraction, and saves the results to CSV files.
%
% INPUT:
% input_wfdb_address - Path to the WFDB file to process
%                      Format: string
%                      Example: 'path/to/record.dat' or 'path/to/record.mat'
%                      Note: The function handles both .dat and .mat WFDB files
%
% ecg_csv_file_name - Path where the ECG features CSV file should be saved
%                     Format: string
%                     Example: 'path/to/output/ecg_features.csv'
%
% lead_names_target - Target lead names to extract from the WFDB file
%                     Format: cell array of strings
%                     Default:  {'I', 'II', 'III', 'aVR', 'aVL', 'aVF', 'V1', 'V2', 'V3', 'V4', 'V5', 'V6'}
%                     Example: {'I', 'II', 'III'}
%
% window_dur - Duration of each processing window in seconds
%              Format: double scalar
%              Default: 60 seconds
%              Range: positive number
%              Note: If 0, processes the entire signal as one window
%
% start_time - Start time for processing in seconds
%              Format: double scalar
%              Default: 0
%              Range: non-negative number
%              Note: If 0, starts from the beginning of the signal
%
% stop_time - Stop time for processing in seconds
%             Format: double scalar
%             Default: 0
%             Range: non-negative number
%             Note: If 0, processes until the end of the signal
%
% flatten_flag - Flag to flatten the ECG features
%                Format:  scalar
%                Default: 0
%                Note: If 1, flattens the all channels ECG features into a single vector
%
% OUTPUT:
% The function saves the following to CSV files:
% 1. ECG features (ecg_csv_file_name)
%    Format: CSV file with columns:
%    - Feature values
%    - Feature names
%    - Feature units
%    - Feature descriptions
%    - Window time information
%    - Lead names

%
% PROCESSING STEPS:
% 1. Load WFDB file and extract signal and metadata
% 2. Preprocess the ECG signal:
%    - Handle NaN values
%    - Apply notch filter (60 Hz)
%    - Remove baseline wander
%    - Apply bandpass filtering
% 3. Extract features in windows:
%    - Process each window separately
%    - Skip windows with NaN values
%    - Extract features using ecg_feature_extraction
% 4. Save results to CSV files
%
% Author: Sajjad Karimi
% Date: Mar 14, 2025
% Location: Emory University, Georgia, USA
% Email: sajjadkarimi91@gmail.com

rng(42);
% SC2I project: waveforms analysis and feature extraction
%
% This function processes WFDB files to extract ECG and PPG features
% and saves them as CSV files.
%
% Inputs:
%   input_wfdb_address - Path to the WFDB file to process
%   save_csv_path - Path where CSV output files should be saved
%
% Dependencies: This script uses some functions from the Open-Source
% Electrophysiological Toolbox (OSET): https://github.com/alphanumericslab/OSET.git

% Load the WFDB file
if strcmp(input_wfdb_address(end-2:end), 'mat')
    [sig, fs, info_ch] = rdmat(input_wfdb_address(1:end-4));
else
    [sig, fs, info_ch] = rddat(input_wfdb_address(1:end-4));
end

if nargin < 3 || isempty(lead_names_target)
    lead_names_target =  {'I', 'II', 'III', 'aVR', 'aVL', 'aVF', 'V1', 'V2', 'V3', 'V4', 'V5', 'V6'};
end

if nargin < 4 || isempty(window_dur)
    window_dur = 60;
end
if nargin < 5 || isempty(start_time)
    start_time = 0;
end

if nargin < 6 || isempty(stop_time)
    stop_time = 0;
end

if nargin < 7 || isempty(flatten_flag)
    flatten_flag = 0;
end

lead_names_temp = cell(1,length(lead_names_target));
ch_I  = zeros(1,length(lead_names_target));
cnt = 0;

% Collect descriptions from info_ch
for ch = 1:length(info_ch)
    if contains(info_ch(ch).Description,lead_names_target)
        cnt = cnt +1;
        ch_I(cnt) = ch;
        lead_names_temp{cnt} = info_ch(ch).Description;
    end
end

if cnt == 0
    error('No ECG channel found in wfdb file')
end

ch_I = ch_I(1:cnt);

lead_names = cell(1,cnt);
for i = 1:cnt
    lead_names{1,i} = lead_names_temp{1,i};
end

ecg_data = sig(:,ch_I)';

% Preprocessing
[C, T] = size(ecg_data); % # of ECG channels
ind_skip = zeros(C, T);

for c = 1:C
    ind_skip(c,:) = isnan(ecg_data(c,:));

    ecg_raw = ecg_data(c,:);
    ecg_raw(isnan(ecg_raw)) = median(ecg_raw,'omitnan');

    % Notch filtering ECG
    fc = 60.0; % powerline frequency
    Qfactor = 45; % Q-factor of the notch filter
    Wo = fc/(fs/2);  BW = Wo/Qfactor; % notch filter parameters
    [b,a] = iirnotch(Wo, BW); % design the notch filter
    ecg_raw = (filtfilt(b, a, ecg_raw(:)))'; % zero-phase non-causal filtering

    % baseline removing
    ecg_raw = ecg_raw - (movmean(movmedian(ecg_raw(:),[round(0.3*fs),round(0.3*fs)]),[round(0.15*fs),round(0.15*fs)]))';

    ecg_raw = ecg_raw - lp_filter_zero_phase(ecg_raw, 0.1/fs); % High pass filter
    ecg_raw = lp_filter_zero_phase(ecg_raw, 30/fs); % Low pass filter

    ecg_data(c,:) = ecg_raw;
end

ind_skip = sum(ind_skip,1)>0;

if stop_time>0 && start_time>0 && stop_time>start_time
    ecg_data = ecg_data(:, start_time*fs+1:stop_time*fs);
    ind_skip = ind_skip(start_time*fs+1:stop_time*fs);

elseif  stop_time==0 && start_time>0
    ecg_data = ecg_data(:, start_time*fs+1:end);
    ind_skip = ind_skip(start_time*fs+1:end);

elseif  stop_time>0 && start_time==0
    ecg_data = ecg_data(:, 1:stop_time*fs);
    ind_skip = ind_skip(1:stop_time*fs);

end

T = size(ecg_data,2); % # of ECG channels

% Feature extraction
if window_dur>0
    win_len = min(T,window_dur*fs);
else
    win_len = T;
end

num_windows = floor(T/win_len);

% % % % Create output filenames
% % % [~, filename, ~] = fileparts(input_wfdb_address);
% % % ecg_csv_file_name = fullfile(save_csv_path, [filename, '_ecg.csv']);

if flatten_flag == 0

    lead_names_csv = cell(cnt,1);
    for j = 1:cnt
        lead_names_csv{j,1} = lead_names{1,j};
    end
else
    lead_names_csv{1} = 'ECG';
end

% Process each window
for n = 1:num_windows
    start_index = (n-1)*win_len+1;
    stop_index = min(n*win_len, T);

    if any(ind_skip(start_index:stop_index))
        continue;
    end

    % ECG feature extraction
    win_ecg_data = ecg_data(:,start_index:stop_index);

    [ecg_features_vector, ecg_feature_info, ecg_fiducial_position, ~] = ...
        ecg_feature_extraction(win_ecg_data, fs, lead_names, [], [], [], [], flatten_flag);

    if all(isnan(ecg_features_vector(:)))
        continue
    end

    if flatten_flag == 0
        window_time_info = repmat([start_index,stop_index], size(ecg_data,1), 1);
        lead_names_csv = cell(cnt,1);
        for j = 1:cnt
            lead_names_csv{j,1} = lead_names{1,j};
        end
    else
        window_time_info = [start_index,stop_index];

    end
    feature_handle_csv(ecg_csv_file_name, ecg_features_vector, ecg_feature_info.names, ...
        ecg_feature_info.units, ecg_feature_info.description, fs, window_time_info, ...
        lead_names_csv, [],  (n==1) );

end

end




