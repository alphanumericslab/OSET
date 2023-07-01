function [b2, b1] = BaseLine2(x, L1, L2, approach)
% Deprecated: BaseLine2 is deprecated. Use baseline_sliding_window instead.
warning('Deprecated: BaseLine2 is deprecated. Use baseline_sliding_window_twice instead.');
[b2, b1] = baseline_sliding_window_twice(x, L1, L2, approach);
