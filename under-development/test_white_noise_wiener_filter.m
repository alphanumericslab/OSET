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
[x, fs] = audioread('/Users/rsameni/Documents/DataFiles/Acoustic_1D_Doppler_OneDrive/121(Good).wav');
x = x(:)';

% Set noise variance and apply Wiener filter
params.nvar = (2e-4)^2;
[y, h] = white_noise_wiener_filter(x, fs, 'fix', params);

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
