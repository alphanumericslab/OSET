function [lowpassed,a] = KLPFilter(y, fc, varargin)
% KLPFilter has been deprecated. Use lp_filter_kalman instead.
warning('KLPFilter has been deprecated. Use lp_filter_kalman instead.');
[lowpassed,a] = lp_filter_kalman(y, fc, varargin{:});