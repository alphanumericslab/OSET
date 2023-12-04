function polarity = ecg_polarity(ecg, fs)
% ECG_POLARITY - Calculate the polarity of multilead ECG data
%
%   polarity = ecg_polarity(ecg, fs)
%
% This function calculates the polarity of multilead ECG data by estimating
% the skewness of the ECG signal after baseline removal and determining
% whether the skewness is positive or negative. The polarity information
% can be used to distinguish QRS complexes with positive and negative
% deflections, or to identify polarity inverted leads (e.g., due to wrong
% lead placement during acquisition)
%
% Inputs:
%   ecg: Multilead ECG data matrix (leads x samples)
%   fs: Sampling frequency of the ECG data
%
% Outputs:
%   polarity: Polarity information for each lead (1 for positive, 0 for negative)
%
% Revision History:
%   2021: First release
%
% Reza Sameni, 2021-2023
% The Open-Source Electrophysiological Toolbox
% https://github.com/alphanumericslab/OSET

fc = 3.0; % An aggressive cut-off frequency for baseline removal

% Apply a low-pass filter to remove the baseline wandering component
baseline = lp_filter_zero_phase(ecg, fc/fs);

% Calculate the skewness of the signal after baseline removal
skw = skew(ecg - baseline);

% Determine the polarity based on the skewness
polarity = skw >= 0;
