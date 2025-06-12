function [feature_vec, feature_info, mean_beat, median_beat] = ecg_snr_features(data, rpeak_indexes)
% snr_features = ecg_snr_features(data, R_peaks_indexes)
% Matlab function for calculating SNR features from ECG signal.
%
% INPUT:
% data          - ECG signal
% R_peaks_indexes - A vector containing the R-peak indices of the ECG signal (expressed as sample points).
%
% OUTPUT:
% snr_features - A struct containing:
%   1. median SNR
%   2. mean SNR
%
% Dependencies:
% 1. `events_snr` function from the OSET package
%
% Author: Seyedeh Somayyeh Mousavi
% Date: Dec 18, 2024
% Location: Emory University, Georgia, USA
% Email: bmemousavi@gmail.com
% Author: Sajjad Karimi
% Location: Emory University, Georgia, USA
% Email: sajjadkarimi91@gmail.com
% Date: Mar 14, 2025

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
% Define feature info
feature_info.names = {'snr_median', 'snr_mean'};
feature_info.units = {'dB', 'dB'};
feature_info.description = {"Median beat SNR", "Mean beat SNR"};

end
