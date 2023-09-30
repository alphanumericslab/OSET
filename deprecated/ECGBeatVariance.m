function [ECGmean, ECGsd] = ECGBeatVariance(x, peaks, wlen)
% ECGBeatVariance has been deprecated. Use avg_beat_calculator instead.
warning('ECGBeatVariance has been deprecated. Use avg_beat_calculator instead.');
[ECGmean, ECGsd] = avg_beat_calculator(x, peaks, wlen);