% ecg_gen_dtw Test Script
%
% This script demonstrates the usage of the ecg_gen_dtw function to generate
% synthetic ECG signals using the warping transformation from the phase domain.
%
% Revision History:
%   2023: First release.
%
% Reza Sameni, 2023
% The Open-Source Electrophysiological Toolbox
% https://github.com/alphanumericslab/OSET

% Clean up environment
close all
clear
clc

% Load ECG model parameters
load('params_00004', 'params');
params = params{1};

% Define parameters for ECG generation
p.alpha = params.a;
p.b = params.b;
p.theta = params.theta;
p.theta0 = 0;
p.f = 90 / 60; % Average heart rate in Hz (beats per second)
p.f_deviations = 0.1; % Percentage of beat-wise heart rate deviations
p.delta_alpha = 0.1; % Percentage of amplitude deviations
p.delta_theta = 0.1; % Percentage of Gaussian center deviations
p.delta_b = 0.1; % Percentage of Gaussian width deviations
warping_order = 2; % Interpolation order for warping transformation

% Simulation parameters
T = 10.0; % Duration of the ECG signal in seconds
fs = 1000; % Sampling rate in Hz
phase_bins = 200; % Number of bins in the phase domain
N = round(T * fs); % Number of samples

% Generate synthetic ECG and phase using ecg_gen_dtw
[ecg, phi] = ecg_gen_dtw(N, fs, p, phase_bins, warping_order);

% Create time vector
time = (0 : N - 1) / fs;

% Plot the synthetic ECG and phase
figure;
yyaxis left
plot(time, ecg)
ylabel('Amplitude (mV)')
yyaxis right
plot(time, phi)
grid
xlabel('Time (s)')
ylabel('Phase (rad)')
set(gca, 'fontsize', 14)
title('Synthetic ECG')
