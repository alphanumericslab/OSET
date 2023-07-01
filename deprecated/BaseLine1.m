function b = BaseLine1(x, L, approach)
% Deprecated: BaseLine1 is deprecated. Use baseline_sliding_window instead.
    warning('Deprecated: BaseLine1 is deprecated. Use baseline_sliding_window instead.');
    b = baseline_sliding_window(x, L, approach);
end