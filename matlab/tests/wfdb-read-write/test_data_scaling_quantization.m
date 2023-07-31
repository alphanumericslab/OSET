clear
clc
close all

% Load the data file
data_struct = load('../../../../../DataFiles/PhysioNetChallenges2023/Sample_files_new/BIDMC_updated/ICARE_0097_20000102_185034.mat');
% data_struct = load('../../../../../DataFiles/PhysioNetChallenges2023/Sample_files_new/ULB_updated/ICARE_0066_20000101_211200.mat');

% Extract the signal data from the loaded structure
data = data_struct.signal;

% Set parameters for data quantization and analysis
params.plot_mode = 'samples';  % 'samples' or 'errors' or 'noplots'
method = 'tanh_sat';  % 'max_scale' or 'tanh_sat' or 'clip'
params.quant_bits = 16;  % 8/16/32/64 fixed-point representation

% Set additional parameters based on the chosen method
if isequal(method, 'tanh_sat')
    params.k_sigma_sat = 20;
end
if isequal(method, 'clip')
    params.clip_level = 10000; % to be defined per case or from device settings
end

% Perform data quantization and analysis
[data_quant, bias, gains, max_amp, snr, snr_med, num_quant_levs] = DataScalingQuantization(data, method, params);

disp(num_quant_levs)
disp(snr)
disp(snr_med)
