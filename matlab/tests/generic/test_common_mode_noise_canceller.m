% test_common_mode_noise_canceller
% 
% This script demonstrates how to use the common_mode_noise_canceller function
% to remove common-mode noise from multichannel data. It imports data from a CSV
% file, processes it to remove common-mode noise, and then plots the original,
% denoised, and smoothed denoised signals for the first 10 channels.
%
% Revision History:
%   2024: First release
%
%   Reza Sameni, 2024
%   The Open-Source Electrophysiological Toolbox
%   https://github.com/alphanumericslab/OSET

clear
close all
clc

% Load data from a CSV file
data = importdata('/Users/rsameni/Documents/caramba.csv');
data = data.data(:, 2:end)'; % Transpose for channels x samples format
fs = 0.5; % Sampling frequency
fc = 0.001; % Cutoff frequency for filtering
itr = 10; % Number of iterations for optimization
lambda = 10; % Regularization parameter
optim_indexes = round(length(data)/4) : round(3*length(data)/4); % Indices for optimization, empty means all

% Call the common_mode_noise_canceller function with specified parameters
[data_den, data_den_smoothed, s, alpha, bl] = common_mode_noise_canceller(data, fs, fc, itr, lambda, optim_indexes, true);

% Plot results 
for ch = 1 : size(data, 1)
    lgnd = {};
    fid = figure('Units', 'inches', 'Position', [0, 0, 18, 6], 'Visible', false); % Set figure size wider
    plot(data(ch, :)); lgnd = cat(1, lgnd, 'Original');
    hold on;
    plot(data_den(ch, :) + bl(ch, :)); lgnd = cat(1, lgnd, 'Denoised');
    plot(data_den_smoothed(ch, :) + bl(ch, :), 'linewidth', 3); lgnd = cat(1, lgnd, 'Smoothed Denoised');
    grid on;
    legend(lgnd, 'interpreter', 'none');
    title(sprintf('Channel %d: Original, Denoised, Smoothed Denoised', ch));
    xlabel('Time (s)');
    ylabel('Amplitude');
    xlim([1, size(data, 2)]);
    set(gca, 'fontsize', 14);
    % Save the figure
    saveas(gcf, sprintf('./results/Channel_%d_Original_Denoised_Smoothed.png', ch), 'png');
    close(fid);
    disp(ch)
end
