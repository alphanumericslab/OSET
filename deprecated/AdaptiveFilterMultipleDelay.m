function [e, n] = AdaptiveFilterMultipleDelay(x, mindelay, maxdelay, taps, mu)
% AdaptiveFilterMultipleDelay has been deprecated. Use adaptive_filter_multi_lag instead.
warning('AdaptiveFilterMultipleDelay has been deprecated. Use adaptive_filter_multi_lag instead.');
[e, n] = adaptive_filter_multi_lag(x, mindelay, maxdelay, taps, mu);