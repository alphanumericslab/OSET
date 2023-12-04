% This file is a demo for testing the Connectivity analysis functions
% This function is produced by Esmaeil Seraj (esmaeil.seraj09@gmail.com)
% 
% Dependencies: Functions "bdf2mat_main.mat", "emg_onset.mat", 
%               "trigger_avg_erp.mat", "BPFilter5.mat", "BaseLine2.mat"
%               and "sig_trend.mat" (provided by the same author). Also 
%               "eeg_read_bdf.mat" (provided by Gleb Tcheslavski
%               (gleb@vt.edu))
% 
% ***NOTE***: Confidential Content, please DO NOT modify or redistribute
%             without permision of the producer. Property of GT-Bionics Lab
%             ECE, Georgia Tech, Atlanta, Georgia, U.S.
%

close all
clear
clc

%% Initialization & conditioning raw signals 
trl = 50;                        % number of trials to be used
emg_ch_ht = 261;                % EMG channel of interest: 261 & 33 & 257
emg_ch_t = 257;                 % EMG channel of interest: 261 & 33 & 257
eeg_ch = 1:16;                  % EEG channel of interest:   6 $ 4 & 4
flag = 'drift';                 % available options are: 'drift', 'nodrift'
filenameformat_ht = 'ht_%d.bdf';
filenameformat_h = 'h_%d.bdf';
filenameformat_t = 't_%d.bdf';

tic
[eeg_data_ht, emg_data_ht, fs, emg_onset_sampl_ht, emg_onset_time_ht] = bdf2mat_main(trl, elec_num, emg_ch_ht, eeg_ch, filenameformat_ht, flag, []);
[eeg_data_h, emg_data_h, ~, emg_onset_sampl_h, emg_onset_time_h] = bdf2mat_main(trl, elec_num, emg_ch_ht, eeg_ch, filenameformat_h, flag);
[eeg_data_t, emg_data_t, ~, emg_onset_sampl_t, emg_onset_time_t] = bdf2mat_main(trl, elec_num, emg_ch_t, eeg_ch, filenameformat_t, flag);
toc

%% Initializing Parameters
pertnum = 20;         % change this to at least 30 iterations in case of having a strong processor 
pairofint = 'all';   % [a, b] form double vector where 'a' and 'b' are pairs of channel for TCPLV estimation
                      % C1: (7) >> C3: (6) >> F3: (3) [change to 'all' for the colored map case]
plot_flag = 'plot';   % decide to visualize the results or not. available options are: 'plot', 'noplot'
freq_rng = [13, 31];  % [a, b] form double vector where 'a' and 'b' are edges of the frequency band of interest
duration = [-3, 2];   % [-a, b] form double vector where 'a' is time required duration prior to movement onset...
                      % and 'b' is the required time duration after the movement onset in seconds

%% Time-Course PLV Dynamics
tic
tcplv_ht = TCPLV(eeg_data_ht, fs, freq_rng, pairofint, emg_onset_time_ht, duration, pertnum);
tcplv_h = TCPLV(eeg_data_h, fs, freq_rng, pairofint, emg_onset_time_h, duration, pertnum);
tcplv_t = TCPLV(eeg_data_t, fs, freq_rng, pairofint, emg_onset_time_t, duration, pertnum);
toc

%% MSC Connectivity Index
tic
% pwcoher = PWCoherence(eeg_data_ht, fs, freq_rng, emg_onset_time_ht, duration, plot_flag);
toc

%% PLV Connectivity Index
tic
% pwplv = PWPLV(eeg_data_ht, fs, emg_onset_time_ht, freq_rng, duration, pertnum, plot_flag);
toc

%% Visualization
% creating scalp-maps
tcplvht_scalp_map = cell(1, size(tcplv_ht, 2));
for i=1:size(tcplv_ht, 2)
    min_data_ht = min(tcplv_ht(:, i));
    interpolate_data_ht = (tcplv_ht(1, i)+tcplv_ht(2, i))/2;
    tcplvht_scalp_map{i} = [min_data_ht, tcplv_ht(1, i), interpolate_data_ht, tcplv_ht(2, i), min_data_ht;...
        min_data_ht, tcplv_ht(3, i), tcplv_ht(4, i), tcplv_ht(5, i), min_data_ht;...
        tcplv_ht(6, i), tcplv_ht(7, i), tcplv_ht(8, i), tcplv_ht(9, i), tcplv_ht(10, i);...
        min_data_ht, tcplv_ht(11, i), tcplv_ht(12, i), tcplv_ht(13, i), min_data_ht;...
        min_data_ht, tcplv_ht(14, i), tcplv_ht(15, i), tcplv_ht(16, i), min_data_ht];
end
tcplvh_scalp_map = cell(1, size(tcplv_h, 2));
for i=1:size(tcplv_h, 2)
    min_data_h = min(tcplv_h(:, i));
    interpolate_data_h = (tcplv_h(1, i)+tcplv_h(2, i))/2;
    tcplvh_scalp_map{i} = [min_data_h, tcplv_h(1, i), interpolate_data_h, tcplv_h(2, i), min_data_h;...
        min_data_h, tcplv_h(3, i), tcplv_h(4, i), tcplv_h(5, i), min_data_h;...
        tcplv_h(6, i), tcplv_h(7, i), tcplv_h(8, i), tcplv_h(9, i), tcplv_h(10, i);...
        min_data_h, tcplv_h(11, i), tcplv_h(12, i), tcplv_h(13, i), min_data_h;...
        min_data_h, tcplv_h(14, i), tcplv_h(15, i), tcplv_h(16, i), min_data_h];
end
tcplvt_scalp_map = cell(1, size(tcplv_t, 2));
for i=1:size(tcplv_t, 2)
    min_data_t = min(tcplv_t(:, i));
    interpolate_data_t = (tcplv_t(1, i)+tcplv_t(2, i))/2;
    tcplvt_scalp_map{i} = [min_data_t, tcplv_t(1, i), interpolate_data_t, tcplv_t(2, i), min_data_t;...
        min_data_t, tcplv_t(3, i), tcplv_t(4, i), tcplv_t(5, i), min_data_t;...
        tcplv_t(6, i), tcplv_ht(7, i), tcplv_ht(8, i), tcplv_ht(9, i), tcplv_ht(10, i);...
        min_data_t, tcplv_t(11, i), tcplv_t(12, i), tcplv_t(13, i), min_data_t;...
        min_data_t, tcplv_t(14, i), tcplv_t(15, i), tcplv_t(16, i), min_data_t];
end

% plotting maps
plotdata_maps = {tcplvht_scalp_map{1}, tcplvht_scalp_map{2}, tcplvht_scalp_map{3}, tcplvht_scalp_map{4}, tcplvht_scalp_map{5},...
    tcplvh_scalp_map{1}, tcplvh_scalp_map{2}, tcplvh_scalp_map{3}, tcplvh_scalp_map{4}, tcplvh_scalp_map{5},...
    tcplvt_scalp_map{1}, tcplvt_scalp_map{2}, tcplvt_scalp_map{3}, tcplvt_scalp_map{4}, tcplvt_scalp_map{5}};
figure
for i=1:length(plotdata_maps)
    subplot(3, 5, i)
    contourf(plotdata_maps{i}, 13)
    colormap jet    
end

% figure
% plot(tcplv_ht, 'b-o','LineWidth', 2)
% hold on; grid on; axis tight;
% plot(tcplv_h, 'r--*','LineWidth', 2)
% plot(tcplv_t, 'g:d','LineWidth', 2)
% legend('hand-tongue', 'hand only', 'tongue only', 'Location', 'northwest')
% xlabel('Time (Sec)'); ylabel('C3-F3 PLV Time-Course Dynamics [8-12Hz]');

%}