function [ECG_mean, ECG_median, ECG_cov, ECG_mean_robust, ECG_median_robust, stacked_beats, beats_sqi] = ecg_features_time_domain(data, peaks, params)
% ecg_features_time_domain - Extracts time-domain ECG mean and median beats, and beat covariance matrix.
%
% Usage:
%   [ECG_mean, ECG_median, ECG_cov, ECG_mean_robust, ECG_median_robust, stacked_beats, beats_sqi] = ecg_features_time_domain(data, peaks, params)
%
% Inputs:
%   data: Matrix of input signals (rows represent channels, and columns represent samples)
%   peaks: Vector containing the indices of detected peaks in the input signal
%   params: Structure with parameters for the time-domain ECG features extraction
%           - params.width: Time window length for feature extraction (optional)
%           - params.window: Window function for covariance and correlation calculations (optional)
%                            If not provided, a rectangular window is used by default.
%
% Outputs:
%   ECG_mean: Time-domain mean ECG beat for each channel (columns)
%   ECG_median: Time-domain median ECG beat for each channel (columns)
%   ECG_cov: Beat covariance matrix for each channel (3rd dimension) in 3D
%   ECG_mean_robust: Time-domain robust mean ECG beat for each channel (columns)
%   ECG_median_robust: Time-domain robust median ECG beat for each channel (columns)
%   stacked_beats: A 3D array of stacked ECG beats for each channel (3rd dimension)
%   beats_sqi: Signal Quality Index (SQI) of each beat for each channel (beats times channels)
% 
% Revision History:
%   2022: First release
%   2023: Renamed from deprecated version TimeDomainECGFeatures
%
% Reza Sameni, 2022-2023 The Open-Source Electrophysiological Toolbox
% https://github.com/alphanumericslab/OSET

I_peaks = find(peaks);
n_beats = length(I_peaks);

% Check if the width parameter is provided, otherwise use the max of RR-intervals in samples by default
if ~isfield(params, 'width') || isempty(params.width)
    params.width = max(diff(I_peaks)); % use max of RR-intervals in samples by default
    if mod(params.width, 2) == 0
        params.width = params.width + 1; % Make window-length odd if even
    end
else
    % Make window-length odd if otherwise
    if mod(params.width, 2) == 0
        params.width = params.width + 1;
        warning('event_width must be odd valued. Automatically modified to the closest greater odd value.');
    end
end

% Check if the window parameter is provided, otherwise use a rectangular window by default
if ~isfield(params, 'window') || isequal(params.window, false)
    window = ones(1, params.width); % rectangular window
else
    window = hamming(params.width)'; % use the specified window function
end

ECG_mean = zeros(params.width, size(data, 1));
ECG_median = zeros(params.width, size(data, 1));
ECG_mean_robust = zeros(params.width, size(data, 1));
ECG_median_robust = zeros(params.width, size(data, 1));
ECG_cov = zeros(params.width, params.width, size(data, 1));
beats_sqi = zeros(n_beats, size(data, 1));
stacked_beats = zeros(n_beats, params.width, size(data, 1));

% Loop through each channel to extract time-domain features
for ch = 1 : size(data, 1)
    stacked_beats(:, :, ch) = event_stacker(data(ch, :), I_peaks, params.width);

    ECG_mean(:, ch) = mean(stacked_beats(:, :, ch), 1); % Compute mean beat
    ECG_median(:, ch) = median(stacked_beats(:, :, ch), 1); % Compute median beat
    ECG_cov(:, :, ch) = cov(stacked_beats(:, :, ch) .* (ones(n_beats, 1) * window)); % Compute beat covariance matrix
    
    % Compute robust mean and median beats using robust_weighted_average function
    [ECG_mean_robust(:, ch), ~, ECG_median_robust(:, ch), ~] = robust_weighted_average(stacked_beats(:, :, ch));

    % Determine the averaging method based on the specified parameter
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
            error('unknown averaging method');
    end

    % Calculate the Signal Quality Index (SQI) for each beat
    for kk = 1 : n_beats
        r = corrcoef(stacked_beats(kk, :, ch) .* window, ECG_avg .* window);
        beats_sqi(kk, ch) = abs(r(1, 2));
    end
end

end
