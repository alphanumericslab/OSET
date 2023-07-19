function [T0, T1] = CalculateTimeLags(peaks, phase)
% CalculateTimeLags has been deprecated. Use calculate_time_lags instead.
warning('CalculateTimeLags has been deprecated. Use calculate_time_lags instead.');
[T0, T1] = calculate_time_lags(peaks, phase);