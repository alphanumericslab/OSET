function [peaks, peak_indexes] = peak_det_simple(x, ff, flag, varargin)
%
% peak_det_simple - A basic R-peak detector
%
%   Syntax: [peaks, peak_indexes] = peak_det_simple(x, ff, flag, omit_close_peaks)
%
%   Inputs:
%       x: Vector of input data
%       ff: Approximate ECG beat-rate in Hertz
%       flag: Search for positive (flag=1), negative (flag=0) peaks, or
%           automatically using the polarity according to the dominant
%           peak/nadir median absolute amplitude
%       omit_close_peaks: omit close peaks after main peak detection
%           (true/1) or not(false/0). Default is 0
%
%   Outputs:
%       peaks: Vector of R-peak impulse train
%       peak_indexes: Vector of R-peak indexes
%
%   Revision History:
%       2006: First release
%       2023: Replaced deprecated version PeakDetection
%       2023: Added mode = 2 to automatically detect sign
%       2023: Added "&& sum(abs(x(index))) > 0" to avoid trivial all-zeros
%             peaks
%
%   Reza Sameni, Davood Fattahi, 2006-2023
%   The Open-Source Electrophysiological Toolbox
%   https://github.com/alphanumericslab/OSET

% Set default values for optional arguments
if nargin > 3 && ~isempty(varargin{1})
    omit_close_peaks = varargin{1};
else
    omit_close_peaks = 0;
end

N = length(x);

rng = floor(0.5 / ff);

peaks = zeros(1, N);
for j = 1 : N
    % Determine the index range for peak search
    if j > rng && j < N - rng
        index = j - rng : j + rng;
    elseif j > rng
        index = N - 2 * rng : N;
    else
        index = 1 : 2 * rng;
    end

    if max(x(index)) == x(j) && sum(abs(x(index))) > 0
        peaks(j) = 1;
    end
end

nadirs = zeros(1, N);
for j = 1 : N
    % Determine the index range for peak search
    if j > rng && j < N - rng
        index = j - rng : j + rng;
    elseif j > rng
        index = N - 2 * rng : N;
    else
        index = 1 : 2 * rng;
    end

    if min(x(index)) == x(j) && sum(abs(x(index))) > 0
        nadirs(j) = 1;
    end
end

switch flag
    case 1
        % No action needed, peaks has already been assigned
    case 0
        % Negative peaks has been explicitly requested
        peaks = nadirs; 
    case 2
        % Automatic mode, select the polarity according to the dominant
        % peak/nadir median absolute amplitude
        if abs(median(x(logical(peaks)))) < abs(median(x(logical(nadirs))))
            peaks = nadirs;
        end
end

% Remove fake peaks
if omit_close_peaks
    I = find(peaks);
    d = diff(I);
    peaks(I(d < rng)) = 0;
end
peak_indexes = find(peaks);
