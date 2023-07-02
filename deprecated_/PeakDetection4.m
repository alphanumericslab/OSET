function [peaks, r] = PeakDetection4(ref, fs, h, fmax, varargin)
% PeakDetection4 has been deprecated. Use peak_detection_matched_filter_robust instead.
warning('PeakDetection4 has been deprecated. Use peak_detection_matched_filter_robust instead.');
[peaks, r] = peak_detection_matched_filter_robust(ref, fs, h, fmax, varargin);


