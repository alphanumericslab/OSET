function test_multi_lead_ecg_viewers
% A script for testing different ECG viewer functions
% 
% Usage:
% testMultiLeadECGPlotter
%
% Revision History:
%   2023: First release
%
% Reza Sameni, 2023
% The Open-Source Electrophysiological Toolbox
% https://github.com/alphanumericslab/OSET

close all
clear
clc

% Define the data file path
datafilepath = '../../../../../DataFiles/Physionet.org/files/ptbdb/1.0.0/';
% directory_list = dir([datafilepath 'patient*']);

% If the data are in WFDB .dat format, convert all .dat files to .mat in bash using this command:
%   find ./*/ -type f -execdir wfdb2mat -r {} \;

% List all .mat files in the specified directory
filelist = dir(fullfile([datafilepath, '**/*.mat']));

% Sampling frequency and mains frequency (Hz)
fs = 1000.0;
f_mains = 60.0;

% Loop through each file
for k = 1 : 1%length(filelist)
    % Load ECG data from the file
    datafilename = [filelist(k).folder '/' filelist(k).name];
    data = load(datafilename);
    data = data.val / 1000.0;

    % Load channel names from the header file
    headerfilename = [filelist(k).folder '/' filelist(k).name];
    header_ID = fopen([headerfilename(1:end-3) 'hea'], 'r');
    line_ = split(fgetl(header_ID), ' ');
    num_ch = str2double(line_{2});
    ch_names = cell(1, num_ch);
    for ch = 1 : num_ch
        line_ = split(fgetl(header_ID), ' ');
        ch_names{ch} = line_{end};
    end
    
    % Limit data to a certain time range
    data = data(:, 1 : round(25.0*fs));

    % Apply bandpass filtering
    data = data - lp_filter_zero_phase(data, 0.25/fs);
    data = lp_filter_zero_phase(data, 80.0/fs);

    % Define parameters for ECG strip viewers
    ref_ch = 1;
    t1_small = 0.0;
    t2_small = 2.5;
    t1_long = 0.0;
    t2_long = 12.0;
    
    % Plot ECG using different viewer functions
    ecg_strip_viewer_multichannel(data(1 : 12, :), fs, ch_names, ref_ch, t1_small, t2_small, t1_long, t2_long, 'Sample Record');
    ecg_strip_viewer_standard_grid(data(1 : 12, 1 : round(t2_long * fs)), fs, ch_names, ref_ch, t1_small, t2_small, 'Sample Record')
    ecg_strip_viewer_waterfall(data(1, :), fs, f_mains, 'Sample Record');
    plot_multichannel_data(data, 4, 'b', fs, 'Sample Record')
end
