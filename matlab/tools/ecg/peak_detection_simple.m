function [peaks, peak_indexes] = peak_detection_simple(x, ff, flag, varargin)
%
% peak_detection_simple - A basic R-peak detector
%
%   Syntax: [peaks, peak_indexes] = peak_detection_simple(x, ff, flag, omit_close_peaks)
%
%   Inputs:
%       x: Vector of input data
%       ff: Approximate ECG beat-rate in Hertz
%       flag: Search for positive (flag=1) or negative (flag=0) peaks
%       omit_close_peaks: omit close peaks after main peak detection
%       (true/1) or not(false/0). Default is 0
%
%   Outputs:
%       peaks: Vector of R-peak impulse train
%       peak_indexes: Vector of R-peak indexes
% 
%   Revision History:
%       2006: First release
%       2023: Extracted from deprecated version PeakDetection()
%
%   Reza Sameni, 2006-2023
%   The Open-Source Electrophysiological Toolbox
%   https://github.com/alphanumericslab/OSET

% Set default values for optional arguments
if nargin > 3 && ~isempty(varargin{1})
    omit_close_peaks = varargin{1};
else
    omit_close_peaks = 0;
end

N = length(x);
peaks = zeros(1, N);

rng = floor(0.5 / ff);

if flag
    for j = 1 : N
        % Determine the index range for peak search
        if j > rng && j < N - rng
            index = j - rng : j + rng;
        elseif j > rng
            index = N - 2 * rng : N;
        else
            index = 1 : 2 * rng;
        end

        if max(x(index)) == x(j)
            peaks(j) = 1;
        end
    end
else
    for j = 1 : N
        % Determine the index range for peak search
        if j > rng && j < N - rng
            index = j - rng : j + rng;
        elseif j > rng
            index = N - 2 * rng : N;
        else
            index = 1 : 2 * rng;
        end

        if min(x(index)) == x(j)
            peaks(j) = 1;
        end
    end
end

% Remove fake peaks
if omit_close_peaks
    I = find(peaks);
    d = diff(I);
    peaks(I(d < rng)) = 0;
end
peak_indexes = find(peaks);
