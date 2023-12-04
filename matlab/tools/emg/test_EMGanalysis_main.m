% This file is a demo for testing the EMG analysis/quantification functions
% This function is produced by Esmaeil Seraj (esmaeil.seraj09@gmail.com)
% 
% Dependencies: Functions "bdf2mat_main.mat", "emg_onset.mat", 
%               "emg_quantification.mat", "BPFilter5.mat", "BaseLine2.mat"
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

%% Initializing and Data Preconditioning
trl = 50;                       % number of trials to be used
emg_ch = 261;                   % EMG channel of interest: 261 & 33 & 257
eeg_ch = 6;                     % EEG channel of interest:   6 $ 4 & 4
flag = 'drift';                 % available options are: 'drift', 'nodrift'
filenameformat_ht = 'ht_%d.bdf';
filenameformat_h = 'h_%d.bdf';
duration = 2.3;

[eeg_data_ht, emg_data_ht, fs, emg_onset_sampl_ht, emg_onset_time_ht] = bdf2mat_main(trl, elec_num, emg_ch_ht, eeg_ch, filenameformat_ht, flag, []);
[eeg_data_h, emg_data_h, ~, emg_onset_sampl_h, emg_onset_time_h] = bdf2mat_main(trl, elec_num, emg_ch_ht, eeg_ch, filenameformat_h, flag);

%% EMG Quantification with LPF+MF Procedure for ECG Cancellation
[emg_avg3_ht, synch_emg2_ht, ecg_estimate2_bl_ht, time_vec_ht] = emg_quantification(emg_data_ht, fs, emg_onset_sampl_ht, duration);
[emg_avg3_h, synch_emg2_h, ecg_estimate2_bl_h, time_vec_h] = emg_quantification(emg_data_h, fs, emg_onset_sampl_h, duration);

%% Measuring the maximum and average onset-period slopes
max_slope_ht = max(diff(emg_avg3_ht));
max_slope_h = max(diff(emg_avg3_h));

onset_delay = [3, 3.1];
onset_vec_ht = emg_avg3_ht(round(onset_delay(1)*fs):round(onset_delay(2)*fs));
avg_slope_ht = mean(diff(onset_vec_ht));
onset_vec_h = emg_avg3_h(round(onset_delay(1)*fs):round(onset_delay(2)*fs));
avg_slope_h = mean(diff(onset_vec_h));

%% visualization
figure
plot((1:length(emg_data_ht{8})-(fs/4))/fs, emg_data_ht{8}(1:end-(fs/4)), 'r', 'LineWidth', 2)
hold on; grid on; axis tight;
plot((1:length(emg_data_ht{8})-(fs/4))/fs, ecg_estimate2_bl_ht{8}(1:end-(fs/4)), 'b', 'LineWidth', 1)
xlabel('Time (Sec)'); ylabel('Amplitude');
legend('Raw EMG Contaminated with ECG',...
    'Estimated ECG through LPF+MF Smoothing', 'location', 'northwest')

figure
plot(time_vec_ht, emg_avg3_ht, 'LineWidth', 1.5)
grid on; axis tight
ylabel(sprintf('EMG Quantity over averaged %d trials HT Task(LPF+MF)', trl)); xlabel('Time (Sec)');

figure
plot(time_vec_h, emg_avg3_h, 'LineWidth', 1.5)
grid on; axis tight
ylabel(sprintf('EMG Quantity over averaged %d trials H Task(LPF+MF)', trl)); xlabel('Time (Sec)');

bar_data = [max_slope_ht max_slope_h; avg_slope_ht avg_slope_h];
figure
bar(bar_data); grid on;
legend('Hand-tongue', 'Hand Only', 'Location', 'Southoutside')
ylabel('First-order time derivatives of quatified EMG curves')
set(gca,'XTick', 1:size(bar_data1, 1))
set(gca,'XTickLabel', {'Maximum Slope During Perfroming the Task',...
    'Average Slope During Muscle Activation [3-3.1 sec]'})
