function [ECG_mean, ECG_median, ECG_cov, ECG_mean_robust, ECG_median_robust, stacked_beats, beats_sqi] = TimeDomainECGFeatures(data, peaks, params)
% Extracting time-domain ECG mean and median beats, and beat covariance matrix.
%
% Usage:
%   [ECG_mean, ECG_median, ECG_cov, ECG_mean_robust, ECG_median_robust, stacked_beats, beats_sqi] = TimeDomainECGFeatures(data, peaks, params)
%
% Reza Sameni
% Oct 2022
% The Open-Source Electrophysiological Toolbox
% https://github.com/alphanumericslab/OSET
%


I_peaks = find(peaks);
n_beats = length(I_peaks);
if ~isfield(params, 'width') || isempty(params.width) % time window length
    params.width = max(diff(I_peaks)); % use max of RR-intervals in samples by default
    if(mod(params.width, 2) == 0)
        params.width = params.width + 1;
    end
else
    % Make window-length odd if otherwise
    if(mod(params.width, 2) == 0)
        params.width = params.width + 1;
        warning('event_width must be odd valued. Automatically modified to the closest greater odd value.');
    end
end

if ~isfield(params, 'window')  || isequal(params.window, false) % window for covariance and correlation calculations
    window = ones(1, params.width);
else
    window = hamming(params.width)';
end
    
ECG_mean = zeros(params.width, size(data, 1));
ECG_median = zeros(params.width, size(data, 1));
ECG_mean_robust = zeros(params.width, size(data, 1));
ECG_median_robust = zeros(params.width, size(data, 1));
ECG_cov = zeros(params.width, params.width, size(data, 1));
beats_sqi = zeros(n_beats, size(data, 1));
for ch = 1 : size(data, 1)
    stacked_beats = EventStacker(data(ch, :), I_peaks, params.width);

    ECG_mean(:, ch) = mean(stacked_beats, 1);
    ECG_median(:, ch) = median(stacked_beats, 1);
    ECG_cov(:, :, ch) = cov(stacked_beats .* (ones(n_beats, 1) * window));
    [ECG_mean_robust(:, ch), ~, ECG_median_robust(:, ch), ~] = RWAverage(stacked_beats);

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

    for kk = 1 : n_beats
        r = corrcoef(stacked_beats(kk, :) .* window, ECG_avg .* window);
        beats_sqi(kk, ch) = abs(r(1, 2));
    end
end

