% A test script for peak_det_acoustic to detect heart beats from a sample
% seismocardiogram record.
%
% This script demonstrates the use of the peak_det_acoustic function to
% process a seismocardiogram (SCG) signal for heart rate analysis. It
% includes loading a sample SCG signal, designing a custom bandpass filter
% to preprocess the signal, setting parameters for peak detection, and
% executing the peak detection process.
%
% Revision History:
%     2024: First release
%
% Reza Sameni, 2024
% The Open-Source Electrophysiological Toolbox
% https://github.com/alphanumericslab/OSET

clear; close all; clc;

% Load a sample SCG signal and its sampling frequency from OSET sample data
load scg_sample.mat samples fs

% Select demo
demo = 'design-filter'; % 'predesigned-filter' or 'design-filter'

switch demo
    case 'predesigned-filter'
        % Option 1: Load pre-designed bandpass filters (uncomment to use)
        load h_bp_filter_4_8hz.mat h_bp

    case 'design-filter'
        % Option 2: Design a custom bandpass filter for SCG signal preprocessing
        % Define filter specifications. All frequency values are in Hz.
        % The following design has been exported from a custom filter design by
        % the MATLAB filterDesigner GUI
        Fstop1 = 4.0;              % First Stopband Frequency
        Fpass1 = 5.0;              % First Passband Frequency
        Fpass2 = 7.0;              % Second Passband Frequency
        Fstop2 = 8.0;              % Second Stopband Frequency
        Dstop1 = 0.01;             % First Stopband Attenuation
        Dpass  = 0.0057563991496;  % Passband Ripple
        Dstop2 = 0.01;             % Second Stopband Attenuation
        dens   = 20;               % Density Factor

        % Calculate the order from the parameters using FIRPMORD.
        [N, Fo, Ao, W] = firpmord([Fstop1 Fpass1 Fpass2 Fstop2]/(fs/2), [0 1 0], [Dstop1 Dpass Dstop2]);

        % Calculate the coefficients using the FIRPM function.
        h_bp = firpm(N, Fo, Ao, W, {dens});
end

% Set parameters for peak_det_acoustic function
params.plot_time_unit = 'seconds';     % Set time unit for plotting
params.h_bp_b = h_bp;                  % Bandpass filter numerator coefficients
params.h_bp_a = 1;                     % Bandpass filter denominator coefficient (FIR filter)
params.filter_type = 'CUSTOM';         % Indicate custom filtering

% Run the peak detection on the SCG sample using the specified parameters
[peak_indexes, hr, hr_smoothed, samples_sat, samples_sat_bp, samples_sat_bp_env] = peak_det_acoustic(samples, fs, params);

% The outputs include:
% peak_indexes - the indices of detected peaks in the signal,
% hr - the raw heart rate calculated from the detected peaks,
% hr_smoothed - the smoothed heart rate signal,
% samples_sat - the signal after amplitude saturation,
% samples_sat_bp - the signal after bandpass filtering,
% samples_sat_bp_env - the power envelope of the bandpass filtered signal.
