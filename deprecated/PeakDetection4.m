function [peaks, r] = PeakDetection4(ref, fs, h, fmax, varargin)
% Deprecated: PeakDetection4 is deprecated. Use peak_detection_matched_filter_robust instead.
warning('Deprecated: PeakDetection4 is deprecated. Use peak_detection_matched_filter_robust instead.');
[peaks, r] = peak_detection_matched_filter_robust(ref, fs, h, fmax, varargin);


