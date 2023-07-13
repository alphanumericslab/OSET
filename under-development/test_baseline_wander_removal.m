% test baseline wander removal by curve fitting
% Status: UNDER TEST
% Reza Sameni, 2005

clear;
close all;
% load LTE07102005_40s.txt -mat;data = LOUETTE_40s';clear LOUETTE_40s
% load fortict_20-Jan-2006.datfortict_20-Jan-2006.dat_acq.mat;
% load fortict5_20-Jan-2006.datfortict5_20-Jan-2006.dat_acq.mat;data = data';
%load Bonhomme3_20-Jan-2006.datBonhomme3_20-Jan-2006.dat_acq.mat;
load patient165_s0323lre; data = data(1:40000, 2)';
% load FOETAL_ECG.dat; data = FOETAL_ECG(:,2:end)';clear FOETAL_ECG;

fs = 1000;

baseline1 = baseline_sliding_window(data, round(fs * 1.0), 'mn');
baseline2 = baseline_sliding_window(data, round(fs * 1.0), 'md');

bl = baseline_sliding_window(data, round(fs * 0.8), 'md');
baseline3 = baseline_sliding_window(bl, round(fs * 1.2), 'mn');
t = (0 : length(data)-1)/fs;

lgnd = {};
figure
plot(t, data); lgnd = cat(1, lgnd, 'raw signal');
hold on
plot(t, baseline1); lgnd = cat(1, lgnd, 'baseline1 mn');
plot(t, baseline2); lgnd = cat(1, lgnd, 'baseline2 md');
plot(t, baseline3); lgnd = cat(1, lgnd, 'baseline3 md-mn');
grid
legend(lgnd)

