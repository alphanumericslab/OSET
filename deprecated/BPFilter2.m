function y = BPFilter2(x, fl, fu, N)
% BPFilter2 has been deprecated. Use bp_filter_dft_varlen instead.
warning('BPFilter2 has been deprecated. Use bp_filter_dft_varlen instead.');
y = bp_filter_dft_varlen(x, fl, fu, N);