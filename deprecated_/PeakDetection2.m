function [peaks, peak_indexes] = PeakDetection2(x, fs, varargin)
% PeakDetection2 has been deprecated. Use peak_detection_modified_pan_tompkins instead.
    warning('PeakDetection has been deprecated. Use peak_detection_modified_pan_tompkins instead.');
    [peaks, peak_indexes] = peak_detection_modified_pan_tompkins(x, fs, varargin{:});
end