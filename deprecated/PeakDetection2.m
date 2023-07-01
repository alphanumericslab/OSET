function [peaks, peak_indexes] = PeakDetection2(x, fs, varargin)
% Deprecated: PeakDetection2 is deprecated. Use peak_detection_modified_pan_tompkins instead.
    warning('Deprecated: PeakDetection is deprecated. Use peak_detection_modified_pan_tompkins instead.');
    [peaks, peak_indexes] = peak_detection_modified_pan_tompkins(x, fs, varargin{:});
end