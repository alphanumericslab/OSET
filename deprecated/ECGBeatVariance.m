function [ECGmean, ECGsd] = ECGBeatVariance(x, peaks, wlen)
% ECGBeatVariance has been deprecated. Use average_beat_calculator instead.
warning('ECGBeatVariance has been deprecated. Use average_beat_calculator instead.');
[ECGmean, ECGsd] = average_beat_calculator(x, peaks, wlen);