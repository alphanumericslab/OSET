function [peaks, peak_indexes] = PeakDetectionAdaptiveHR(x, ff, fs, varargin)
% PeakDetectionAdaptiveHR has been deprecated. Use peak_detection_adaptive_hr instead.
warning('PeakDetectionAdaptiveHR has been deprecated. Use peak_detection_adaptive_hr instead.');
[peaks, peak_indexes] = peak_detection_adaptive_hr(x, ff, fs, varargin{:});
