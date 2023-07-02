function [y1, y2, Pbar, Phat, PSmoothed, Kgain] = KalmanARFilter(x, b, a, varargin)
% KalmanARFilter has been deprecated. Use kalman_ar_filter instead.
warning('KalmanARFilter has been deprecated. Use kalman_ar_filter instead.');
[y1, y2, Pbar, Phat, PSmoothed, Kgain] = kalman_ar_filter(x, b, a, varargin);