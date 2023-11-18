clear; close all; clc

fs_input = 11025;
fs_output = 4000.0;
T = 0.1;
N = round(T * fs_input);
t_input = (0 : N - 1)/fs_input;

% Sample test of the Resampler class using the mexFunction interface
signal_input = sin(2 * pi * 157.13 * t_input) + sin(2 * pi * 256.3 * t_input);
conversion_rate = fs_input/fs_output; % Resample to half the input sample rate (input/output sampling frequency)
order = 4; % Use fourth-order interpolation

% Call the resampler_mex function to perform resampling
signal_resampled = resampler_mex(signal_input, conversion_rate, order);
t_output = (0 : length(signal_resampled) - 1)/fs_output;

% Call the resampler_mex function again to perform resampling back to the
% original sampling frequency
signal_input_reconstructed = resampler_mex(signal_resampled, 1./conversion_rate, order);
t_input_reconstructed = (0 : length(signal_input_reconstructed) - 1)/fs_input;

% Plot the original and resampled signals
figure
plot(t_input, signal_input, 'b-', 'linewidth', 2);
hold on
plot(t_output, signal_resampled, 'r--', 'linewidth', 2);
plot(t_input_reconstructed, signal_input_reconstructed, 'g:', 'linewidth', 2);
legend('Input Signal', 'Resampled Signal', 'Reconstructed Signal');
grid

input_vs_reconstructed_delay = finddelay(signal_input , signal_input_reconstructed)