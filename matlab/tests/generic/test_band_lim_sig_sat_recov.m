% A test script for band_lim_sig_sat_recov(), which recovers saturated signals by iterative band-limited filtering
% Reza Sameni, 2024
% The Open-Source Electrophysiological Toolbox (OSET)
% https://github.com/alphanumericslab/OSET

close all
clear
clc

% Load signal from file (to be replaced with an arbitray signal) 
datafilename = 'input_audio_file_name';
[x, fs] = audioread(datafilename, 'native');
x = double(x(:)'); % Convert to double

n_itr = 10; % Number of iterations
th = 28000; % Saturation threshold
fu = 750.0 / fs; % High cutoff frequency normalized by sampling frequency
fl = [] / fs; % High cutoff frequency normalized by sampling frequency
th_bottom_end = []; % Lower bound for saturated values after recovery. By defult: `th_bottom - 0.25 * abs(th_bottom)`
th_top_end = []; % Upper bound for saturated values after recovery. By default: `th_top + 0.25 * abs(th_top)`
verbose = true; 
y = band_lim_sig_sat_recov(x, -th, th, n_itr, fl, fu, [], [], verbose);

% PLOT
t = (0 : length(x) - 1) / fs;
figure
plot(t, y, 'DisplayName', 'Filtered signal')
hold on
plot(t, x, 'DisplayName', 'Input signal')
grid on
legend show
xlabel('Time (s)')
ylabel('Amplitude')
set(gca, 'fontsize', 16)
