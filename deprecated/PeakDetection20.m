function [peaks, peak_indexes] = PeakDetection20(x, ff, th, varargin)
% PeakDetection20 has been deprecated. Use peak_det_amp_threshold instead.
warning('PeakDetection20 has been deprecated. Use peak_det_amp_threshold instead.');
if nargin > 3 && ~isempty(varargin{1})
    flag = varargin{1};
else
    flag = 2;
end
level = 'MEDIAN';
[peaks, peak_indexes] = peak_det_amp_threshold(x, ff, th, flag, level);
