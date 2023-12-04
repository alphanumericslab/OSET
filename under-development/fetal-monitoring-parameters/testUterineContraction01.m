% A test file for Uterine Contraction
%
% Reza Sameni
% Aug 2020

clear;
clc;
close all;

% fs = 4000; % Sampling frequency
% alpha = -0.1;
% UC_iden_wlen = 5.0;
% UC_avg_wlen = 15.0;
% Tachysystolic_th = 40;
% EHG_fl = 2.0;
% EHG_fu = 400;
% EHG_env_win_param = 0.02; % envelope calculation window length
% data = load("/Users/rsameni/Documents/GitRepository/NasimKatebiThesis/DUS/SampleDoppler/03008294-d6c4-3ac7-ad3f-de4d43b6bfa7Patient-2.mat");
% x = data.wavRecordings{1}';
% x = x(1:round(40.0*fs));

% fs = 200; % sampling frequency
f0 = 50; % powerline frequency
fmat = 1.0; % maternal heart rate in Hz
alpha = 0.15;
UC_iden_wlen = 50.0;
UC_avg_wlen = 150.0;
Tachysystolic_th = 40;
EHG_fl = 1.0;
EHG_fu = 50;
EHG_env_win_param = 0.05; % envelope calculation window length


% 1- Load data:
% data = load('ice001_l_1of1m.mat');
% data = data.val(:, round(5*fs):end-round(5*fs));

% data = load('ice002_p_1of3m.mat');
% data = data.val(:, round(5*fs):round(2500*fs)-1);

fs = 20; % sampling frequency
% data = load('tpehgt_t007m.mat');
data = load('tpehgt_p006m.mat');
data = data.val;

data_len = size(data, 2);
% 2- Remove baseline:
data_HP = data - LPFilter(data, 1.0/fs);

% figure
% hold on
% plot(data_HP(1, :));
% plot(data(7, :));
% grid
wwlen = round(30.0 * fs);
figure
subplot(211);
plot((0:data_len-1)/fs/60., data(7, :));
hold on
plot((0:data_len-1)/fs/60., data_HP(1, :));
grid
subplot(212);
% spectrogram(data_HP(1, :), hamming(wwlen), round(wwlen/2), 512, fs, 'yaxis');
spectrogram(data_HP(1, :), hamming(wwlen), wwlen - 2, 512, fs, 'yaxis');


% 3- Multi-stage BSS to remove maternal ECG:

% ICA using JADE algorithm
W = jadeR(data_HP);
indp_sources = W * data_HP;

% Find maternal R-peaks and apply PiCA
chmat = 1;
mpeaks = PeakDetection(indp_sources(chmat, :), 1.0/fs);
I_peaks = find(mpeaks);
[pica_sources, W, A] = PiCA(data_HP, mpeaks);

% Plot the maternal R-peaks to make sure the right component was chosen and processed
figure
hold on
plot(indp_sources(chmat, :))
plot(I_peaks, indp_sources(chmat, I_peaks), 'ro')
grid
title('The first maternal ECG source and its R-peaks');

% Remove the maternal ECG:
% all_channels_mean = RWAverage(data_HP);
% data_BP_denoised = data_BP - all_channels_mean(ones(1, size(data_BP, 1)), :);

% data_BP_denoised = A(:, 3:end) * s(3:end, :);
% s(1 : 2, :) = 0;
pica_sources_den = pica_sources;
pica_sources_den(1, :) = pica_sources(1, :) - wden(pica_sources(1, :), 'rigrsure', 's','mln', 2, 'coif5');
pica_sources_den(2, :) = pica_sources(2, :) - wden(pica_sources(2, :), 'rigrsure', 's','mln', 2, 'coif5');
data_HP_denoised = A * pica_sources_den;

% Plot the first twi periodic components to make sure the right component was chosen and processed
figure
subplot(211);
plot(pica_sources(1, :));
hold on
plot(pica_sources_den(1, :));
grid
subplot(212);
plot(pica_sources(2, :));
hold on
plot(pica_sources_den(2, :));
grid

% Wo = f0/(fs/2);
% BW = Wo/35;
% [num, den] = iirnotch(Wo, BW); 
% data_HP_denoised_notched = filtfilt(num, den, data_HP_denoised);

% data_BP = LPFilter(data_HP_notched, 30.0/fs);
% data_preprocessed = data_HP_denoised_notched;
% data_preprocessed = data_HP_denoised;
data_preprocessed = zeros(size(data_HP_denoised));
for kk = 1 : size(data_HP_denoised, 1)
%     data_preprocessed(kk, :) = bandpass(data_HP_denoised(kk, :), [35, 65], fs);
    data_preprocessed(kk, :) = wden(data_HP_denoised(kk, :), 'rigrsure', 's','sln', 6, 'coif5');
    data_preprocessed(kk, :) = bandpass(data_preprocessed(kk, :), [20, 80], fs);    
end



ch_spectral = 1;
wwlen = round(1.0 * fs);
figure
subplot(211);
plot((0:data_len-1)/fs, data_HP_denoised(ch_spectral, :));
hold on
plot((0:data_len-1)/fs, data_preprocessed(ch_spectral, :));
grid
subplot(212);
spectrogram(data_preprocessed(ch_spectral, :), hamming(wwlen), round(wwlen/2) - 2, 512, fs, 'yaxis');

% PlotECG(data_BP_denoised, 4, 'b', fs);
% 
% figure
% hold on
% plot(data_HP');
% plot(all_channels_mean, 'k', 'linewidth', 2);
% grid;

% y = data - LPFilter(data, 1.0/fs);
% y = data_BP_denoised - LPFilter(data_BP_denoised, 1.0/fs);

% W = jadeR(y);
% s = W * y;
% x = sqrt(sum(s(3, :).^2, 1));
% xx = s(3, :);

% xx = sqrt(sum(data_preprocessed.^2, 1));
% xx = abs(data_preprocessed(1, :));

UA_channel = 1;
wlen = round(EHG_env_win_param * fs); % envelope detection window length in samples
x = sqrt(filtfilt(ones(1, wlen), wlen, data_preprocessed(UA_channel, :).^2)); % IUP signal

figure
hold on
plot((0:data_len-1)/fs, data_preprocessed(UA_channel, :));
plot((0:data_len-1)/fs, x);
grid

T = length(x);
t = (0 : T - 1)/fs;
plotflag = 1;
% [UC, UC_onsets, UC_counts, UC_counts_avg, Tachysystolic] = UterineContractions(x, 'EHG', fs, alpha, UC_iden_wlen, UC_avg_wlen, Tachysystolic_th, EHG_fl, EHG_fu, EHG_env_win_param);
[UC, UC_onsets, UC_counts, UC_counts_avg, Tachysystolic] = UterineContractions(x, 'IUP', fs, alpha, UC_iden_wlen, UC_avg_wlen, Tachysystolic_th, plotflag);

% figure
% hold on
plot(t, UC);
stem(t(UC_onsets), UC_onsets(UC_onsets));
plot(t, UC_counts);
plot(t, UC_counts_avg);
plot(t, Tachysystolic, 'linewidth', 2);
% grid