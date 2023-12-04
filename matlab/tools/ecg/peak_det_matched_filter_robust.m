function [peaks, r] = peak_det_matched_filter_robust(ref, fs, h, fmax, varargin)
% peak_det_matched_filter_robust: R-peak detector based on a matched filter
%
% Syntax: [peaks, r] = peak_det_matched_filter_robust(ref, fs, h, fmax, itr)
%
% Inputs:
%   ref:    Vector of input data.
%   fs:     Sampling rate.
%   h:      Template waveform.
%   fmax:   Maximum expected frequency of the R-peaks.
%   itr:    Number of iterations to run the post-matched filter peak detector (Default: 1, run peak detector only once)
%
% Outputs:
%   peaks:  Vector of R-peak impulse train.
%   r:      Filtered output after matched filtering.
%
%   Revision History:
%       2008: First release
%       2023: Renamed from deprecated version PeakDetection4()
%
%   Reza Sameni, 2006-2023
%   The Open-Source Electrophysiological Toolbox
%   https://github.com/alphanumericslab/OSET

% Parse the vararg input
if isempty(varargin)
    itr = 1; % Default value if itr
else
    itr = varargin{1};
end

N = length(ref);
L = length(h);

h = h(end:-1:1); % Reverse the template waveform

w = floor(L/2);

r = filter(h, 1, [ref zeros(1, w-1)]); % Matched filtering
r = r(w : N+w-1); % Trim the filtered output to match the length of ref

% Use peak_det_local_search function to detect peaks
peaks = peak_det_local_search(r, fmax/fs, 1, itr);
