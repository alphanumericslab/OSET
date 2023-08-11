function [peaks, r] = PeakDetection4(ref, fs, h, fmax, varargin)
% PeakDetection4 has been deprecated. Use peak_det_matched_filter_robust instead.
warning('PeakDetection4 has been deprecated. Use peak_det_matched_filter_robust instead.');
[peaks, r] = peak_det_matched_filter_robust(ref, fs, h, fmax, varargin);


