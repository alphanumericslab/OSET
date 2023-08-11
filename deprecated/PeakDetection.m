function [peaks, peak_indexes] = PeakDetection(x, ff, varargin)
% PeakDetection has been deprecated. Use peak_det_local_search instead.
warning('PeakDetection has been deprecated. Use peak_det_local_search instead.');
[peaks, peak_indexes] = peak_det_local_search(x, ff, varargin{:});
