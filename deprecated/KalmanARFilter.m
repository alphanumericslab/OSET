function [y1, y2, Pbar, Phat, PSmoothed, Kgain] = KalmanARFilter(x, b, a, varargin)
% KalmanARFilter has been deprecated. Use ar_filter_kalman instead.
warning('KalmanARFilter has been deprecated. Use ar_filter_kalman instead.');
[y1, y2, Pbar, Phat, PSmoothed, Kgain] = ar_filter_kalman(x, b, a, varargin{:});