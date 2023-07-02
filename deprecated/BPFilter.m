function y = BPFilter(x, fl, fu)
% BPFilter has been deprecated. Use bp_filter_fft instead.
warning('BPFilter has been deprecated. Use bp_filter_fft instead.');
y = bp_filter_fft(x, fl, fu);