function y = BPFilter(x, fl, fu)
% BPFilter has been deprecated. Use bp_filter_dft instead.
warning('BPFilter has been deprecated. Use bp_filter_dft instead.');
y = bp_filter_dft(x, fl, fu);