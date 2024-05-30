% A script to test the vcg_construction.m function for generating VCG with
% physical coordinates from multilead ECG data

clear; close all; clc;

datafilepath = '../../../../../DataFiles/physionet.org/files/ptbdb/1.0.0/'; fs = 1000.0;

filelist = dir(fullfile([datafilepath, '**/*.mat']));  % get list of all mat files
fl = 0.25;
fh = 80.0;
% f_hr = 1.4; % default initial heart rate value in Hz
% top_source_num = 1;
% num_rounds = 2; % default number of R-peak detection rounds; see peak_det_local_search() help for details
% baseline_fraction = 0.1; % percentage length of the beat (from the start) considered as the baseline segment (roughly)
% avg_beat_right_expansion_frac = 1.2;
% nsca_method = 'COV'; % 'COV' or 'CORR'
% refch = 1;
% data_augmented = [];

for k = 73 : 73%length(filelist)
    datafilename = [filelist(k).folder '/' filelist(k).name];
    data = load(datafilename);
    data = data.val;
    headerfilename = [filelist(k).folder '/' filelist(k).name];
    header_ID = fopen([headerfilename(1:end-3) 'hea'], 'r');
    line_ = split(fgetl(header_ID), ' ');
    num_ch = str2double(line_{2});
    ch_names = cell(1, num_ch);
    for ch = 1 : num_ch
        line_ = split(fgetl(header_ID), ' ');
        ch_names{ch} = line_{end};
    end

    data = data - lp_filter_zero_phase(data, fl/fs);
    data = lp_filter_zero_phase(data, fh/fs);

    ecg = data(1:12, :);
    vcg = data(13:15, :);

    vcg_rec = vcg_construction(ecg, ch_names, [], true);

    t = (0 : size(data, 2)-1)/fs;  
    for ch = 1 : 3
        figure
        hold on
        plot(t, vcg(ch, :))
        plot(t, vcg_rec(ch, :))
        grid
        legend('Actual Frank leads', 'Reconstructed VCG')
        xlabel('time[s]')
        ylabel('Amplitude[mV]')
        title(ch_names{ch + 12})
        set(gca, 'fontsize', 16)
    end

end
