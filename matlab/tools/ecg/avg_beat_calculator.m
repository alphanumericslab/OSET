function [ECGmean, ECGsd] = avg_beat_calculator(x, peaks, wlen, varargin)
% avg_beat_calculator - Calculation of the mean and standard deviation of ECG waveforms in different beats.
%
% Usage:
%   [ECGmean, ECGsd] = avg_beat_calculator(x, peaks, wlen, offset_method)
%
% Inputs:
%   x: Input ECG signal
%   peaks: ECG peaks
%   wlen: Window length around the peaks (from each side of the peak)
%   offset_method: Offset method for the beats before averaging
%       'none': No offsetting
%       'mean': Offset by the mean value of the beat
%       'peak': Offset by the peak value of the beat
%
% Outputs:
%   ECGmean: Mean ECG beat
%   ECGsd: Standard deviation of ECG beats
%
% Revision History:
%   2010: First release
%   2023: Renamed from deprecated version ECGBeatVariance and integrated
%       with offset feature from depricated function ECGBeatVariance2
%
% Reza Sameni, 2010-2023
% The Open-Source Electrophysiological Toolbox
% https://github.com/alphanumericslab/OSET

% Check for offset method argument
if nargin > 3 && ~isempty(varargin{1})
    offset_method = varargin{1};
    if ~isequal(offset_method, 'none') && ~isequal(offset_method, 'mean') && ~isequal(offset_method, 'peak')
        error('Undefined offset method. Please use ''none'', ''mean'', or ''peak''.');
    end
else
    offset_method = 'none';
end

% Add padding to the input signal and peaks
x = [zeros(1, wlen), x, zeros(1, wlen)];
peaks = [zeros(1, wlen), peaks, zeros(1, wlen)];

% Find the indices of peaks
I = find(peaks);

N = length(I);

ECG = zeros(N, 2 * wlen);
for i = 1:N
    switch offset_method
        case 'none'
            offset = 0;
        case 'mean'
            offset = mean(x(I(i) - wlen + 1:I(i) + wlen));
        case 'peak'
            offset = x(I(i));
    end
    ECG(i, :) = x(I(i) - wlen + 1:I(i) + wlen) - offset;
end

ECGmean = mean(ECG, 1);
ECGsd = std(ECG, 1, 1);
end
