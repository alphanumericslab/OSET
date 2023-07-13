function y = BPFilter6(x, fc, bw, factor)
% BPFilter6 has been deprecated. Use bp_filter_complex_equiripple instead.
warning('BPFilter6 has been deprecated. Use bp_filter_complex_equiripple instead.');
y = bp_filter_complex_equiripple(x, fc, bw, factor);