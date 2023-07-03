function y = BPFilter4(x,fc,bw,order)
% BPFilter4 has been deprecated. Use bp_filter_complex_ma_alt instead.
warning('BPFilter4 has been deprecated. Use bp_filter_complex_ma_alt instead.');
y = bp_filter_complex_ma_alt(x,fc,bw,order);