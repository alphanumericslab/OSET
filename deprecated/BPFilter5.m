function y = BPFilter5(x, fc, bw, order)
% BPFilter5 has been deprecated. Use bp_filter_complex_ma instead.
warning('BPFilter5 has been deprecated. Use bp_filter_complex_ma instead.');
y = bp_filter_complex_ma(x, fc, bw, order);