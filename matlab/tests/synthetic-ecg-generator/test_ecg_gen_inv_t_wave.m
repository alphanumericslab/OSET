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
% load('params_00012', 'params');
% params = params{1};

load('params_01000');
params = params{3};

% Define parameters for ECG generation
p.alpha = params.a;
p.b = params.b;
p.theta = params.theta;
p.theta0 = 0;
p.f = 70 / 60; % Average heart rate in Hz (beats per second)
p.f_deviations = 0.1; % Percentage of beat-wise heart rate deviations
p.delta_alpha = 0.1; % Percentage of amplitude deviations
p.delta_theta = 0.1; % Percentage of Gaussian center deviations
p.delta_b = 0.1; % Percentage of Gaussian width deviations
p.T_inv_shape_factor = 0.1; % T-wave inversion shape factor (determines tanh inversion function slope)
p.T_onset = 0.02; % Anticipated onset of the T-wave 
p.T_offset = inf; % Anticipated offset of the T-wave
warping_order = 2; % Interpolation order for warping transformation
seed = 0; % Random number generator seed

% Simulation parameters
T = 10.0; % Duration of the ECG signal in seconds
fs = 1000; % Sampling rate in Hz
phase_bins = 300; % Number of bins in the phase domain
N = round(T * fs); % Number of samples

% Generate synthetic ECG and phase using ecg_gen_warping
no_inversion = true;
[ecg_ref, ~] = ecg_gen_inv_t_wave(N, fs, p, phase_bins, warping_order, seed, no_inversion);
[ecg, phi] = ecg_gen_inv_t_wave(N, fs, p, phase_bins, warping_order, seed);

% Create time vector
time = (0 : N - 1) / fs;

% Plot the synthetic ECG and phase
figure;
yyaxis left
hold on
plot(time, ecg_ref)
plot(time, ecg)
legend('Reference ECG', 'Inverted T-wave')
ylabel('Amplitude (mV)')
yyaxis right
plot(time, phi)
grid
xlabel('Time (s)')
ylabel('Phase (rad)')
set(gca, 'fontsize', 14)
title('Synthetic ECG')
