function y = BPFilterComplex(x, w0, N)
% BPFilterComplex has been deprecated. Use bp_filter_complex_ma_fixed_tap instead.
warning('BPFilterComplex has been deprecated. Use bp_filter_complex_ma_fixed_tap instead.');
y = bp_filter_complex_ma_fixed_tap(x, w0, N);