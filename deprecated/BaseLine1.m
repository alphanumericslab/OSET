function b = BaseLine1(x, L, approach)
% BaseLine1 has been deprecated. Use baseline_sliding_window instead.
warning('BaseLine1 has been deprecated. Use baseline_sliding_window instead.');
b = baseline_sliding_window(x, L, approach);
