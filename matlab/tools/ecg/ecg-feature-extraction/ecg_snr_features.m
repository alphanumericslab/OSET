function [feature_vec, feature_info, mean_beat, median_beat] = ecg_snr_features(data, rpeak_indexes)

% function [feature_vec, feature_info] = ecg_snr_features(data, R_peaks_indexes)
% Extract features from ECG using beat signal-to-noise ratio (SNR)
%
% Inputs:
%   data: ECG signal as a vector (in microvolts) (1D array).
%   R_peaks_indexes:  A vector containing the R-peak indices of the ECG signal (expressed as sample points).
%
% Outputs:
% feature_vec: A vector contains
%   1. median SNR
%   2. mean SNR
%
%  feature_info: A structure contains feature descriptions, names and
%  units.
%
%  mean_beat: Average ECG beat signal (1D array)
%
% Dependencies:
%   1. `events_snr` function from the OSET package
%
% Author:
%   Seyedeh Somayyeh Mousavi
%   Sajjad Karimi
%   Reza Sameni
%   Emory University, Georgia, USA
%   Email: bmemousavi@gmail.com
%   First: Date: SEP 24, 2024
%   Second: Date: AUG 3, 2025

%% Constants
increasing_beat_length = 1.2;

% Calculate beat length
beat_length = round(increasing_beat_length * median(diff(rpeak_indexes)));
if mod(beat_length, 2) == 0
    beat_length = beat_length + 1;
end
beat_length = ceil([0.3*beat_length,0.7*beat_length]);
% Calculate SNR values
[snr_median, snr_mean, mean_beat, median_beat] = events_snr(data, rpeak_indexes, beat_length);

% Assign the median and mean SNR values to the struct
snr_features.median = median(snr_median);
snr_features.mean = median(snr_mean);

feature_vec = [snr_features.median, snr_features.mean];

% Convert Inf to NaN
feature_vec(isinf(feature_vec)) = NaN;

% Define feature info
feature_info.names = {'snr_median', 'snr_mean'};
feature_info.units = {'dB', 'dB'};
feature_info.description = {"Median beat SNR", "Mean beat SNR"};

end
