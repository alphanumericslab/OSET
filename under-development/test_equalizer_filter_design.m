clear
close all
clc

fs = 1000;
spectral_window = 3.0;
spectral_overlap = 2.75;
spectral_len = round(spectral_window*fs);

ecg1 = load('SampleECG1kHz1.mat');
ecg1 = ecg1.data;
ecg1 = ecg1 - mean(ecg1);
Sx = pwelch(ecg1, hamming(spectral_len), round(spectral_overlap*fs),spectral_len, fs, 'twosided');

ecg2 = load('SampleECG1kHz2.mat');
ecg2 = ecg2.data(:, 2:end);
ecg2 = ecg2' - mean(ecg2)';

Sy = zeros(size(ecg2, 1), spectral_len);
for ch = 1 : size(ecg2, 1)
    Sy(ch, :) = pwelch(ecg2(ch, :), hamming(spectral_len), round(spectral_overlap*fs),spectral_len, fs, 'twosided')';
end

Sy = median(Sy, 1);

params.epsilon                   = eps;
params.fs                        = fs;
params.filter_len                = 101;
params.lambda                    = 10.0;
params.smooth_spectrum           = false;
params.innovation_filter_type    = 'LINEAR_PHASE'; % 'LINEAR_PHASE', 'MIN_PHASE'
params.plot_results              = true;


h = equalizer_filter_design(Sx, Sy, params);
ecg1_equalized = filter(h, sqrt(sum(h.^2)), ecg1);

Sx_equalized = pwelch(ecg1_equalized, hamming(spectral_len), round(spectral_overlap*fs),spectral_len, fs, 'twosided');

ff = fs * (0:spectral_len-1) / spectral_len;
figure
lgnd = {};
plot(ff - fs/2, fftshift(10*log10(abs(Sx))), 'linewidth', 3); hold on; lgnd{end+1} = 'Sx';
plot(ff - fs/2, fftshift(10*log10(abs(Sy))), 'linewidth', 3); hold on; lgnd{end+1} = 'Sy';
plot(ff - fs/2, fftshift(10*log10(abs(Sx_equalized))), 'linewidth', 3); hold on; lgnd{end+1} = 'Sx_equalized';
grid
legend(lgnd, 'interpreter', 'none')

time1 = (0:length(ecg1)-1)/fs;
time2 = time1 - length(h)/2/fs;
figure
plot(time1, ecg1)
hold on
plot(time2, ecg1_equalized)
grid
axis tight