function testMultiLeadECGPlotter

close all
clear
clc

% Convert all dat files in bash using this command: find ./*/ -type f -execdir wfdb2mat -r {} \;

datafilepath = '../../../../../DataFiles/Physionet.org/files/ptbdb/1.0.0/';
% directory_list = dir([datafilepath 'patient*']);
% filelist = dir(fullfile([datafilepath, '**/*lr*.mat']));  % get list of all mat files
filelist = dir(fullfile([datafilepath, '**/*.mat']));  % get list of all mat files
fs = 1000.0;

for k = 1 : 1%length(filelist)
    datafilename = [filelist(k).folder '/' filelist(k).name];
    data = load(datafilename);
    data = data.val / 1000.0;
    
    headerfilename = [filelist(k).folder '/' filelist(k).name];
    header_ID = fopen([headerfilename(1:end-3) 'hea'], 'r');
    line_ = split(fgetl(header_ID), ' ');
    num_ch = str2double(line_{2});
    ch_names = cell(1, num_ch);
    for ch = 1 : num_ch
        line_ = split(fgetl(header_ID), ' ');
        ch_names{ch} = line_{end};
    end

    data = data(:, 1 : round(25.0*fs));
    
    % bandpass filter
    data = data - lp_filter_zero_phase(data, 0.25/fs);
    data = lp_filter_zero_phase(data, 80.0/fs);
    
    % plot the multichannel ECG
    ref_ch = 1;
    t1_small = 0.0;
    t2_small = 3.0;
    t1_long = 0.0;
    t2_long = 10.0;
    MultiLeadECGPlotter(data(1 : 12, :), ch_names, fs, ref_ch, t1_small, t2_small, t1_long, t2_long, 'Record #0000');
    % MultiLeadECGPlotterStandardSize(data(1 : 12, 1 : round(t2_long * fs)), ch_names, fs, ref_ch, t1_small, t2_small, 'Record #0000');
    
end
