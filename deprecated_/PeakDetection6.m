function [peaks, peak_indexes] = PeakDetection6(x,ff,th,varargin)
% PeakDetection6 has been deprecated. Use peak_detection_amp_threshold instead.
    warning('PeakDetection6 has been deprecated. Use peak_detection_local_search instead.');
[peaks, peak_indexes] = peak_detection_amp_threshold(x, ff, th, varargin{:});