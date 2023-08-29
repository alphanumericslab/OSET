% ecg_gen_twa Test Script
%
% This script demonstrates the usage of the ecg_gen_twa function to generate
% synthetic ECG signals with T-wave alternans (TWA) using the warping
% transformation from the phase domain.
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

% Load pre-defined ECG model parameters (for demonstration purposes)
load('params_01000');
params = params{3};

% Define parameters for ECG generation
p.alpha = params.a;
p.b = params.b;
p.theta = params.theta;
p.theta0 = 0;
p.f = 70 / 60; % Average heart rate in Hz (beats per second)
p.f_deviations = 0.1; % Percentage of beat-wise heart rate deviations
p.delta_alpha = 0.05; % Percentage of amplitude deviations
p.delta_theta = 0.05; % Percentage of Gaussian center deviations
p.delta_b = 0.05; % Percentage of Gaussian width deviations
p.TWA_shape_factor = 0.1; % T-wave segment shape factor (determines tanh inversion function slope)
p.T_onset = 0.02; % Anticipated onset of the T-wave 
p.T_offset = inf; % Anticipated offset of the T-wave
p.TWA_fluct_factor = 0.1; % Fluctuation factor for TWA
warping_order = 2; % Interpolation order for warping transformation
seed = 0; % Random number generator seed
p.STM = [0 1; 1 0]; % Exact alternation: the T-wave will alternate in each beat
% p.STM = [.2 .8; .8 .2]; % Probabilistic alternation with high probability of state transition
% p.STM = [.9 .1; .99 .01]; % Abnormalities with low probability of occurrence
p.state0 = 1; % First T-wave type (1 or 2)

% Simulation parameters
T = 15.0; % Duration of the ECG signal in seconds
fs = 1000; % Sampling rate in Hz
phase_bins = 300; % Number of bins in the phase domain
N = round(T * fs); % Number of samples

% Generate synthetic ECG and phase using ecg_gen_twa
[ecg, phi] = ecg_gen_twa(N, fs, p, phase_bins, warping_order, seed);

% Create time vector
time = (0 : N - 1) / fs;

% Plot the synthetic ECG and phase
figure;
hold on
plot(time, ecg)
ylabel('Amplitude (mV)')
xlabel('Time (s)')
grid
set(gca, 'fontsize', 14)
title('Synthetic ECG with T-wave alternans')
