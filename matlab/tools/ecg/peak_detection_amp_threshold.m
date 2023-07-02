function [peaks, peak_indexes] = peak_detection_amp_threshold(x, ff, th, varargin)
% peak_detection_amp_threshold - R-peak detector based on max search and level thresholding
%
% Syntax: [peaks, peak_indexes] = peak_detection_amp_threshold(x, ff, th, varargin)
%
% Inputs:
%   x: Vector of input data.
%   ff: Approximate ECG beat-rate in Hertz, normalized by the sampling frequency.
%   th: Peaks smaller than this fraction of the max peak amplitude are neglected.
%   flag (optional): Search for positive (flag=1) or negative (flag=0) peaks. By default,
%               the maximum absolute value of the signal determines the peak sign.
%
% Outputs:
%   peaks: Vector of R-peak impulse train.
%   peak_indexes: Indexes of the detected R-peaks.
%
% Notes:
% - The R-peaks are found from a peak search in windows of length N, where
%   N corresponds to the R-peak period calculated from the given ff. R-peaks
%   with periods smaller than N/2 or greater than N are not detected.
% - The signal baseline wander is recommended to be removed before the
%   R-peak detection.
%
%   Revision History:
%       2006: First release
%       2020: Used logical indexing for speed improvement
%       2023: Renamed from deprecated version PeakDetection6()
%
%   Reza Sameni, 2006-2023
%   The Open-Source Electrophysiological Toolbox
%   https://github.com/alphanumericslab/OSET

if nargin == 4
    flag = varargin{1};
else
    flag = abs(max(x)) > abs(min(x));
end

omit_close_peaks = 0;
[peaks, ~] = peak_detection_simple(x, ff, flag, omit_close_peaks);

% Remove small peaks below the threshold
I = find(peaks);
mx = max(abs(x(I)));
J = abs(x(I)) < th * mx;
peaks(I(J)) = 0;

peak_indexes = find(peaks);
