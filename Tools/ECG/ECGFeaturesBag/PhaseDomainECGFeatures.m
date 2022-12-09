function [ECG_mean, ECG_median, ECG_cov, ECG_mean_robust, ECG_median_robust, stacked_beats, beats_sqi] = PhaseDomainECGFeatures(data, peaks, params)
% Extracting phase-domain ECG mean and median beats, and beat covariance matrix.
%
% Usage:
%   [ECG_mean, ECG_median, ECG_cov, ECG_mean_robust, ECG_median_robust, stacked_beats, beats_sqi] = PhaseDomainECGFeatures(data, peaks, params)
% 
% Reza Sameni
% Oct 2022
% The Open-Source Electrophysiological Toolbox
% https://github.com/alphanumericslab/OSET
%


I_peaks = find(peaks);
if ~isfield(params, 'bins') || isempty(params.bins) % number of phase bins
    params.bins = median(diff(I_peaks)); % use median of RR-intervals in samples by default
end

n_full_beats = length(I_peaks) - 1;

if ~isfield(params, 'window') || isempty(params.window) || isequal(params.window, false)% window for covariance and correlation calculations
    window = ones(1, params.bins);
else
    window = hamming(params.bins)';
end

ECG_mean = zeros(params.bins, size(data, 1));
ECG_median = zeros(params.bins, size(data, 1));
ECG_mean_robust = zeros(params.bins, size(data, 1));
ECG_median_robust = zeros(params.bins, size(data, 1));
ECG_cov = zeros(params.bins, params.bins, size(data, 1));
beats_sqi = zeros(n_full_beats, size(data, 1));
for ch = 1 : size(data, 1)
    stacked_beats = zeros(n_full_beats, params.bins);
    for k = 2 : n_full_beats + 1
        previous_peak = I_peaks(k - 1);
        current_peak = I_peaks(k);
        beat = data(ch, previous_peak : current_peak);
        x_phase_warped = LinearWarp(beat, params.bins + 1);
        stacked_beats(k - 1, :) = x_phase_warped(1 : end - 1);
    end

    stacked_beats = fftshift(stacked_beats, 2);
    ECG_mean(:, ch) = mean(stacked_beats, 1);
    ECG_median(:, ch) = median(stacked_beats, 1);
    ECG_cov(:, :, ch) = cov(stacked_beats .* (ones(n_full_beats, 1) * window));
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

    for kk = 1 : n_full_beats
        r = corrcoef(stacked_beats(kk, :), ECG_avg);
        beats_sqi(kk, ch) = abs(r(1, 2));
    end

end