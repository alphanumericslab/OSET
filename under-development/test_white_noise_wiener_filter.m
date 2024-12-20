% Test script for white_noise_wiener_filter (white noise suppression using
% Wiener filter)
% 
% Reza Sameni, 2024
% The Open-Source Electrophysiological Toolbox
% https://github.com/alphanumericslab/OSET

close all
clc
clear

% Load the input audio signal
[x, fs] = audioread('/Users/rsameni/Documents/DataFiles/Acoustic_1D_Doppler_OneDrive/104(Good).wav', 'native');
% [x, fs] = audioread('/Users/rsameni/Documents/DataFiles/Acoustic_1D_Doppler_OneDrive/46(Good).wav', 'native');
% [x, fs] = audioread('/Users/rsameni/Documents/DataFiles/Acoustic_1D_Doppler_OneDrive/121(Good).wav', 'native');
% [x, fs] = audioread('/Users/rsameni/Documents/DataFiles/Acoustic_1D_Doppler_OneDrive/127(Intermediate).wav', 'native');
% [x, fs] = audioread('/Users/rsameni/Documents/DataFiles/Acoustic_1D_Doppler_OneDrive/128(Poor).wav', 'native');

x = double(x(:)');

% Set noise variance and apply Wiener filter
params.f_cut = 1500.0;
[y, h] = white_noise_wiener_filter(x, fs, 'fixed-freq', params);
% params.nvar = (5)^2;
% [y, h] = white_noise_wiener_filter(x, fs, 'fix', params);

% Alternative example using automatic noise variance estimation
% [y, h] = white_noise_wiener_filter(x, fs, 'min-spectral-power');

% Calculate and display the Signal-to-Noise Ratio (SNR)
noise = x - y;
snr = 10 * log10(mean(y.^2) / mean(noise.^2));
disp(['SNR = ', num2str(snr)])

% Plot input and filtered signals
t = (0 : length(x)-1) / fs;
figure
plot(t, x, 'DisplayName', 'Input signal')
hold on
plot(t, y, 'DisplayName', 'Filtered signal')
grid on
legend show
xlabel('Time (s)')
ylabel('Amplitude')
set(gca, 'fontsize', 16)
