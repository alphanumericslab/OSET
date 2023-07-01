function [peaks, peak_indexes] = PeakDetection(x, ff, varargin)
% Deprecated: PeakDetection is deprecated. Use peak_detection_local_search instead.
    warning('Deprecated: PeakDetection is deprecated. Use peak_detection_local_search instead.');
    [peaks, peak_indexes] = peak_detection_local_search(x, ff, varargin{:});
end