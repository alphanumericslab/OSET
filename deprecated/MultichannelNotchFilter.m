function [s, W, A, B, C] = MultichannelNotchFilter(x, fc, Q, fs)
% MultichannelNotchFilter has been deprecated. Use sca_notch_filter instead.
warning('MultichannelNotchFilter has been deprecated. Use sca_notch_filter instead.');
[s, W, A, B, C] = sca_notch_filter(x, fc, Q, fs);