% Test script for ecg_gen_warping
% 
% This script demonstrates the usage of the ecg_gen_warping function to generate
% synthetic ECG signals using the warping transformation from the phase to time domain.
% 
% Revision History:
%   2023: First release.
%
% Reza Sameni, 2023
% The Open-Source Electrophysiological Toolbox
% https://github.com/alphanumericslab/OSET

close all
clear
clc

% Load parameters
load('params_00001', 'params'); % Load ECG parameter settings
params = params{1}; % Assume params contains a cell array with a single structure

% Define parameters for ECG generation
p.alpha = params.a;
p.b = params.b;
p.theta = params.theta;
p.f = 75/60; % Average heart rate in Hz (BPM/60)
p.f_deviations = 0.1; % Percentage of beat-wise heart rate deviations (Hz)
p.delta_alpha = 0.1;
p.delta_theta = 0.1;
p.delta_b = 0.1;
warping_order = 3; % Interpolation order for warping transformation

T = 10.0; % Desired duration of the synthetic ECG signal (seconds)
fs = 1000; % Sampling rate of the synthetic ECG signal (Hz)
phase_bins = 200; % Number of bins in the phase domain
N = round(T * fs); % Number of samples

% Generate synthetic ECG signal using ecg_gen_warping
[ecg, phi] = ecg_gen_warping(N, fs, p, phase_bins, warping_order);

% Create a time vector
time = (0 : N-1) / fs;

% Plot the synthetic ECG signal and corresponding phase
figure;
yyaxis left
plot(time, ecg)
ylabel('Amplitude (mV)')
yyaxis right
plot(time, phi)
grid on
xlabel('Time (s)')
ylabel('Phase (rad)')
set(gca, 'fontsize', 14)
title('Synthetic ECG')
legend('ECG', 'Phase')
