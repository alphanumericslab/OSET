% test_linear_phase_equalizer_design.m
% Design and visualize a linear phase equalizer filter
% 
%   Reza Sameni, 2024
%   The Open-Source Electrophysiological Toolbox
%   https://github.com/alphanumericslab/OSET

% Context: let's say we have a signal with spectrum S_source(f), which we
% want to equalize to a target spectrum S_target(f). From stochastic
% process theory, this can be achived with a filter with squared magnitude
% response: |H(f)|^2 = S_target(f)./S_source(f). The following can be
% performed to design the filter's magnitude response |H(f)|:

% Clear the environment
clear; close all; clc;

fs = 48000; % Sampling frequency in Hz

% Step 1: define the desired frequency response (by visualizing |H(f)|)

% Normalized frequency bands (0 to Nyquist) - change these according to |H(f)|
f = [0 100 200 1000 2000 24000]/(fs/2);

% Desired amplitude response in each band - change these according to |H(f)|
% a = [1 1 2 2 0.5 0.5]; % in linear scale
a = 10.^([0 0 6 6 -6 -6] / 20); % in dB scale

% Step 2: design the filter
N = 500; % Filter order (change as per need)
% FIR filter coefficients using frequency sampling
b = fir2(N, f, a); 

[H, W] = freqz(b,1,512,fs); % the filter's frequency response

% Alternatively, using least squares FIR filter design
% b = firls(N, f, a);

% Step 3: apply the filter on sample signals

t = 0:1/fs:1-1/fs; % time vector

% Sine wave:
% x = sin(2*pi*500*t) + 0.5*sin(2*pi*1500*t); % Example signal with two frequencies (change as needed)

% Random input:
x = randn(1, length(t));

y = filter(b, 1, x); % Apply the equalizer

% Visualize the results
figure
subplot(211)
plot(W, 20*log10(abs(H)))
xlim([W(1), W(end)])
grid on
title('Frequency response of the equalizer')
xlabel('Frequency (Hz)')
ylabel('amplitude');
set(gca, 'fontsize', 16)

subplot(212)
[Pxx, Wxx] = pwelch(y, [], [], [], fs);
plot(Wxx, 10*log10(abs(Pxx)))
grid on
title('Sample signal''s spectra')
xlabel('Frequency (Hz)')
ylabel('amplitude');
set(gca, 'fontsize', 16)

figure; 
subplot(2,1,1)
plot(t, x)
title('Original signal')
xlabel('time (s)')
ylabel('amplitude')
grid 
set(gca, 'fontsize', 16)
subplot(2,1,2)
plot(t, y)
title('Filtered signal')
xlabel('time (s)')
ylabel('amplitude')
grid
set(gca, 'fontsize', 16)

