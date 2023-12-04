function y = BPFilter3(x, fc, bw, wname)
% BPFilter3 has been deprecated. Use bp_filter_dft_windowed instead.
warning('BPFilter3 has been deprecated. Use bp_filter_dft_windowed instead.');
y = bp_filter_dft_windowed(x, fc, bw, wname);