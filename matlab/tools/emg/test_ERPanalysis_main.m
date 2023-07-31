% This file is a demo for testing the ERD/ERS analysis functions
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
trl = 50;                       % number of trials to be used
emg_ch_ht = 261;                % EMG channel of interest: 261 & 33 & 257
emg_ch_t = 257;                 % EMG channel of interest: 261 & 33 & 257
eeg_ch = 6;                     % EEG channel of interest:   6 $ 4 & 4
flag = 'drift';                 % available options are: 'drift', 'nodrift'
elec_num = 16;
filenameformat_ht = 'ht_%d.bdf';
filenameformat_h = 'h_%d.bdf';
filenameformat_t = 't_%d.bdf';

[eeg_data_ht, emg_data_ht, fs, emg_onset_sampl_ht, emg_onset_time_ht] = bdf2mat_main(trl, elec_num, emg_ch_ht, eeg_ch, filenameformat_ht, flag, []);
[eeg_data_h, emg_data_h, ~, emg_onset_sampl_h, emg_onset_time_h] = bdf2mat_main(trl, elec_num, emg_ch_ht, eeg_ch, filenameformat_h, flag);
[eeg_data_t, emg_data_t, ~, emg_onset_sampl_t, emg_onset_time_t] = bdf2mat_main(trl, elec_num, emg_ch_t, eeg_ch, filenameformat_t, flag);

%% ERP time-course estimation through triggered-average of EEG ensembles
duration = 2;    % required signal duration after movement onset in seconds
freq_band = 'beta';                  % available options: 'alpha' or 'beta'
ref_per = [-1, -0.3]; % double vector in [-a, -b] form where -a and -b are the edges of reference segment
cof_intv = 3;                             % confidence interval coefficient

[erp_ht, synch_eeg_ht, trigger_time_sec_ht, time_vec_ht, synch_emg_ht] = trigger_avg_erp(eeg_data_ht, fs, freq_band, emg_onset_time_ht, duration, emg_data_ht);
[erp_h, synch_eeg_h, trigger_time_sec_h, time_vec_h, synch_emg_h] = trigger_avg_erp(eeg_data_h, fs, freq_band, emg_onset_time_h, duration, emg_data_h);
[erp_t, synch_eeg_t, trigger_time_sec_t, time_vec_t, synch_emg_t] = trigger_avg_erp(eeg_data_t, fs, freq_band, emg_onset_time_t, duration, emg_data_t);

[erd_area_ht, ers_area_ht, quant_erp_ht] = erp_quantification(erp_ht, fs, trigger_time_sec_ht, ref_per, cof_intv);
[erd_area_h, ers_area_h, quant_erp_h] = erp_quantification(erp_h, fs, trigger_time_sec_h, ref_per, cof_intv);
[erd_area_t, ers_area_t, quant_erp_t] = erp_quantification(erp_t, fs, trigger_time_sec_t, ref_per, cof_intv);

%% ERP Time/Frequency Representation
method = 'STFT';               % available options: 'NBCH', 'STFT' or 'CWT'
% [erp_tf, ~, ~, ~, freq_vec_tf, time_vec_tf] = trigger_avg_TF_erp(eeg_data_ht, emg_data_ht, fs, emg_onset_sampl_ht, duration, method);

%% Visualization
figure
subplot(511)
plot(time_vec_ht, synch_emg_ht{6})
title('Synchronized EMG Signals Based on Trigger Time')
axis([min(time_vec_ht) max(time_vec_ht) -1.5e4 1.5e4])
ylabel('Trl (6)')
grid on
subplot(512)
plot(time_vec_ht, synch_emg_ht{7})
axis([min(time_vec_ht) max(time_vec_ht) -1.5e4 1.5e4])
ylabel('Trl (7)')
grid on
subplot(513)
plot(time_vec_ht, synch_emg_ht{8})
axis([min(time_vec_ht) max(time_vec_ht) -1.5e4 1.5e4])
ylabel('Trl (8)')
grid on
subplot(514)
plot(time_vec_ht, synch_emg_ht{9})
axis([min(time_vec_ht) max(time_vec_ht) -1.5e4 1.5e4])
ylabel('Trl (9)')
grid on
subplot(515)
plot(time_vec_ht, synch_emg_ht{10})
axis([min(time_vec_ht) max(time_vec_ht) -1.5e4 1.5e4])
ylabel('Trl (10)')
xlabel('time (sec)')
grid on

figure
plot(time_vec_ht, quant_erp_ht, 'LineWidth', 1.5)
title('ERP Curve Extracted By Trigger-Averaged EEG Trials')
line([trigger_time_sec_ht, trigger_time_sec_ht], [0, max(quant_erp_ht)],  'Color', [1, 0.4, 0], 'linewidth', 1.5, 'linestyle', '-.')
axis tight; grid on
ylabel('Event Related Potential'); xlabel('time (sec)')

bar_data1 = [erd_area_ht{1} erd_area_h{1} erd_area_t{1};...
    ers_area_ht{1} ers_area_h{1} ers_area_t{1}];
bar_data2 = [erd_area_ht{2} 0 0;...
    ers_area_ht{2} ers_area_h{2} 0];
bar_data3 = [0 0 0;...
    ers_area_ht{3} 0 0];
figure
subplot(131)
bar(bar_data1); grid on;
legend('Hand-tongue', 'Hand Only', 'Tongue Only', 'Location', 'Southoutside')
ylabel('Normalized Area Enclosed by ERP Curve and Mean \pm STD')
set(gca,'XTick', 1:size(bar_data1, 1))
set(gca,'XTickLabel', {'ERD Events #1', 'ERS Events #1'})
subplot(132)
bar(bar_data2); grid on;
legend('Hand-tongue', 'Hand Only', 'Tongue Only', 'Location', 'Southoutside')
ylabel('Normalized Area Enclosed by ERP Curve and Mean \pm STD')
set(gca,'XTick', 1:size(bar_data1, 2))
set(gca,'XTickLabel', {'ERD Events #2', 'ERS Events #2'})
subplot(133)
bar(bar_data3); grid on;
legend('Hand-tongue', 'Hand Only', 'Tongue Only', 'Location', 'Southoutside')
ylabel('Normalized Area Enclosed by ERP Curve and Mean \pm STD')
set(gca,'XTick', 1:size(bar_data1, 3))
set(gca,'XTickLabel', {'ERD Events #3', 'ERS Events #3'})

bar_new_data_erd = [erd_area_ht{1} erd_area_h{1} erd_area_t{1};...
    erd_area_ht{2} 0 0; 0 0 0];
bar_new_data_ers = [ers_area_ht{1} ers_area_h{1} ers_area_t{1};...
    ers_area_ht{2} ers_area_h{2} 0; ers_area_ht{3} 0 0];
figure
subplot(211)
bar(bar_new_data_erd); grid on;
legend('Hand-tongue', 'Hand Only', 'Tongue Only')
ylabel('Normalized ERD Areas')
set(gca,'XTick', 1:size(bar_data1, 1))
set(gca,'XTickLabel', {'ERD Events #1', 'ERD Events #2', 'ERD Events #3'})
subplot(212)
bar(bar_new_data_ers); grid on;
legend('Hand-tongue', 'Hand Only', 'Tongue Only')
ylabel('Normalized ERS Areas')
set(gca,'XTick', 1:size(bar_data1, 2))
set(gca,'XTickLabel', {'ERS Events #1', 'ERS Events #2', 'ERS Events #3'})

tt_ht = (trigger_time_sec_ht*fs)-((1.5)*fs):(trigger_time_sec_ht*fs)+((1.5)*fs);
tt_h = (trigger_time_sec_h*fs)-((1.5)*fs):(trigger_time_sec_h*fs)+((1.5)*fs);
tt_t = (trigger_time_sec_t*fs)-((1.5)*fs):(trigger_time_sec_t*fs)+((1.5)*fs);
figure
subplot(131)
plot(time_vec_ht, quant_erp_ht, 'LineWidth', 1.5)
xlabel('Time (Sec)'); ylabel('Relative ERP Magnitude (%) [Hand-Tongue]')
grid on; axis tight;
subplot(132)
plot(time_vec_h, quant_erp_h, 'LineWidth', 1.5)
xlabel('Time (Sec)'); ylabel('Relative ERP Magnitude (%) [Hand Only]')
grid on; axis tight;
subplot(133)
plot(time_vec_t, quant_erp_t, 'LineWidth', 1.5)
xlabel('Time (Sec)'); ylabel('Relative ERP Magnitude (%) [Tongue Only]')
grid on; axis tight;

