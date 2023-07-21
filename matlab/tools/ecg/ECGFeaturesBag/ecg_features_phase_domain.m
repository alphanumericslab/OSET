function [ECG_mean, ECG_median, ECG_cov, ECG_mean_robust, ECG_median_robust, stacked_beats_all_channels, beats_sqi] = ecg_features_phase_domain(data, peaks, params)
% ecg_features_phase_domain - Extracts phase-domain ECG mean and median beats, and beat covariance matrix.
%
% Usage:
%   [ECG_mean, ECG_median, ECG_cov, ECG_mean_robust, ECG_median_robust, stacked_beats, beats_sqi] = ecg_features_phase_domain(data, peaks, params)
%
% Inputs:
%   data: Matrix of input signals (rows represent channels, and columns represent samples).
%   peaks: Vector containing the indices of detected peaks in the input signal.
%   params: A structure with optional parameters:
%       - bins: Number of phase bins for the phase-domain representation.
%       - window: Window for covariance and correlation calculations (if empty or false, it uses a rectangular window).
%                If specified, it should be a row vector of the same length as the number of bins.
%
% Outputs:
%   ECG_mean: Phase-domain ECG mean beats (rows represent phase bins, and columns represent channels).
%   ECG_median: Phase-domain ECG median beats (rows represent phase bins, and columns represent channels).
%   ECG_cov: Phase-domain ECG beat covariance matrices (size: bins x bins x channels).
%   ECG_mean_robust: Robust phase-domain ECG mean beats (rows represent phase bins, and columns represent channels).
%   ECG_median_robust: Robust phase-domain ECG median beats (rows represent phase bins, and columns represent channels).
%   stacked_beats_all_channels: 3D matrix containing phase-domain stacked beats for all channels (size: n_full_beats x bins x channels).
%   beats_sqi: Squared correlation coefficient between each beat and the corresponding averaged beat (size: n_full_beats x channels).
%
% Reference for the notion of ECG phase domain:
%   R. Sameni, C. Jutten, and M. B. Shamsollahi. Multichannel
%       electrocardiogram decomposition using periodic component analysis. IEEE
%       Transactions on Biomedical Engineering, 55(8):1935-1940, Aug. 2008.
%
% Revision History:
%   2022: First release
%   2023: Renamed from deprecated version PhaseDomainECGFeatures
%
% Reza Sameni, 2022-2023 The Open-Source Electrophysiological Toolbox
% https://github.com/alphanumericslab/OSET

% Get the indices of detected peaks
I_peaks = find(peaks);
if ~isfield(params, 'bins') || isempty(params.bins) % number of phase bins
    params.bins = median(diff(I_peaks)); % use median of RR-intervals in samples by default
end

n_full_beats = length(I_peaks) - 1;

if ~isfield(params, 'window') || isempty(params.window) || isequal(params.window, false) % window for covariance and correlation calculations
    window = ones(1, params.bins);
else
    window = hamming(params.bins)';
end

% Initialize output variables
ECG_mean = zeros(params.bins, size(data, 1));
ECG_median = zeros(params.bins, size(data, 1));
ECG_mean_robust = zeros(params.bins, size(data, 1));
ECG_median_robust = zeros(params.bins, size(data, 1));
ECG_cov = zeros(params.bins, params.bins, size(data, 1));
beats_sqi = zeros(n_full_beats, size(data, 1));
stacked_beats_all_channels = zeros(n_full_beats, params.bins, size(data, 1));

% Loop through each channel
for ch = 1 : size(data, 1)
    % Preallocate matrix for stacked phase-domain beats
    stacked_beats = zeros(n_full_beats, params.bins);
    
    % Iterate through each full beat to compute the phase-domain representation
    for k = 2 : n_full_beats + 1
        previous_peak = I_peaks(k - 1);
        current_peak = I_peaks(k);
        beat = data(ch, previous_peak : current_peak);
        x_phase_warped = LinearWarp(beat, params.bins + 1);
        stacked_beats(k - 1, :) = x_phase_warped(1 : end - 1);
    end

    % Shift the phase-domain representation
    stacked_beats = fftshift(stacked_beats, 2);
    
    % Compute mean, median, covariance, and robust mean/median
    ECG_mean(:, ch) = mean(stacked_beats, 1);
    ECG_median(:, ch) = median(stacked_beats, 1);
    ECG_cov(:, :, ch) = cov(stacked_beats .* (ones(n_full_beats, 1) * window));
    [ECG_mean_robust(:, ch), ~, ECG_median_robust(:, ch), ~] = robust_weighted_average(stacked_beats);

    % Store stacked beats for all channels
    stacked_beats_all_channels(:, :, ch) = stacked_beats;

    % Determine the averaging method based on the given parameter
    switch params.BEAT_AVG_METHOD
        case 'MEAN'
            ECG_avg = ECG_mean(:, ch)';
        case 'MEDIAN'
            ECG_avg = ECG_median(:, ch)';
        case 'ROBUST_MEAN'
            ECG_avg = ECG_mean_robust(:, ch)';
        case 'ROBUST_MEDIAN'
            ECG_avg = ECG_median_robust(:, ch)';
        otherwise
            error('Unknown averaging method');
    end

    % Compute squared correlation coefficient (SQI) between each beat and the corresponding averaged beat
    for kk = 1 : n_full_beats
        r = corrcoef(stacked_beats(kk, :), ECG_avg);
        beats_sqi(kk, ch) = abs(r(1, 2))^2;
    end

end
