function [ECGmean, ECGsd] = ECGBeatVariance2(x, peaks, wlen)
% ECGBeatVariance2 has been deprecated. Use avg_beat_calculator instead.
warning('ECGBeatVariance2 has been deprecated. Use avg_beat_calculator instead.');
offset_method = 'peak';
[ECGmean, ECGsd] = avg_beat_calculator(x, peaks, wlen, offset_method);