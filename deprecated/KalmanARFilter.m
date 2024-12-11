function [y1, y2, Pbar, Phat, PSmoothed, Kgain] = KalmanARFilter(x, b, a, varargin)
% KalmanARFilter has been deprecated. Use arma_kalman_filter instead.
warning('KalmanARFilter has been deprecated. Use arma_kalman_filter instead.');
[y1, y2, ~, Pbar, Phat, PSmoothed, Kgain] = arma_kalman_filter(x, b, a, varargin{:});