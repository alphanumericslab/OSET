function [e,n] = AdaptiveFilter(x, delay, taps, mu)
% AdaptiveFilter has been deprecated. Use adaptive_filter instead.
warning('AdaptiveFilter has been deprecated. Use adaptive_filter instead.');
[e, n] = adaptive_filter(x, delay, taps, mu);