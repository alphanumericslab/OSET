function [peaks, peak_indexes] = peak_det_adaptive_hr(x, ff, fs, varargin)
% peak_det_adaptive_hr - R-peak detector based on max search over sliding window with adaptive width
% Syntax:
%   [peaks, peak_indexes] = peak_det_adaptive_hr(x, ff, fs, flag, [th, th2], rejection_threshold)
% 
% Method:
%   detects R-peaks in the ECG signal x using a max search algorithm over a sliding window
%   with adaptive width. The function returns the binary impulse train of detected peaks
%   (PEAKS) and their corresponding indexes (PEAK_INDEXES) in the input signal x.
%
% Inputs:
%   x: Vector of the input ECG signal.
%   ff: Approximate ECG beat-rate in Hz.
%   fs: Sampling frequency in Hz.
%   flag (optional): Search for positive (flag=1) or negative (flag=0) peaks. 
%       By default, the maximum absolute value of the signal determines the peak sign.
%   [th, th2] (optional): Thresholds used for the first and second round of peak detection (default: [0.5, 0.8]).
%   rejection_threshold (optional): Threshold for rejecting fake peaks (default: 0.3).
%
% Outputs:
%   peaks: Binary impulse train indicating the detected R-peaks.
%   peak_indexes: Indexes of the detected R-peaks in the input signal.
%
% Note:
%   - The signal baseline wander should be removed before R-peak detection.
%
% Revision History:
%   2021: First release.
%   2023: Renamed from deprecated version PeakDetectionAdaptiveHR.
%
% Reza Sameni, 2021-2023
% The Open-Source Electrophysiological Toolbox
% https://github.com/alphanumericslab/OSET

N = length(x);
peaks0 = zeros(1, N);
peaks = zeros(1, N);

% Check optional input arguments
if nargin > 3 && ~isempty(varargin{1})
    flag = varargin{1};
else
    flag = abs(max(x)) > abs(min(x));
end

if nargin > 4 && ~isempty(varargin{2})
    [th, th2] = varargin{2};
else
    th = 0.5;
    th2 = 0.8;
end

if nargin > 5 && ~isempty(varargin{3})
    rejection_threshold = varargin{3};
else
    rejection_threshold = 0.3; % rejects peaks below 25% of the median peak amplitudes
end

rng = floor(th / (ff / fs));

% First round peak detection
if flag
    for j = 1 : N
        if j > rng && j < N - rng
            index = j - rng : j + rng;
        elseif j > rng
            index = N - 2 * rng : N;
        else
            index = 1 : 2 * rng;
        end
        
        if max(x(index)) == x(j)
            peaks0(j) = 1;
        end
    end
else
    for j = 1 : N
        if j > rng && j < N - rng
            index = j - rng : j + rng;
        elseif j > rng
            index = N - 2 * rng : N;
        else
            index = 1 : 2 * rng;
        end
        
        if min(x(index)) == x(j)
            peaks0(j) = 1;
        end
    end
end

% Smooth the RR intervals
II = find(peaks0);
RR_intervals = diff(II);
RR_intervals_smoothed = baseline_sliding_window(RR_intervals, 3, 'md');
RR_intervals_smoothed2 = baseline_sliding_window(RR_intervals_smoothed, 3, 'mn');
ff = fs ./ [RR_intervals_smoothed2(1), RR_intervals_smoothed2];

% Handle edge cases
if II(1) > 1
    II = [1, II];
    ff = [mean(ff(1:2)), ff];
end
if II(end) < N
    II = [II, N];
    ff = [ff, mean(ff(end-1:end))];
end

% Interpolate the beat-rate for all samples
ff_interpolated = interp1(II, ff, 1:N);
rng2 = floor(th2 ./ (ff_interpolated / fs));

% Second round peak detection
if flag
    for j = 1 : N
        if j > rng2(j) && j < N - rng2(j)
            index = j - rng2(j) : j + rng2(j);
        elseif j > rng2(j)
            index = N - 2 * rng2(j) : N;
        else
            index = 1 : 2 * rng2(j);
        end
        
        if max(x(index)) == x(j)
            peaks(j) = 1;
        end
    end
else
    for j = 1 : N
        if j > rng2(j) && j < N - rng2(j)
            index = j - rng2(j) : j + rng2(j);
        elseif j > rng2(j)
            index = N - 2 * rng2(j) : N;
        else
            index = 1 : 2 * rng2(j);
        end
        
        if min(x(index)) == x(j)
            peaks(j) = 1;
        end
    end
end

% Remove fake peaks
I = find(peaks);
peak_amps = median(x(I));
J = abs(x(I)) < rejection_threshold * abs(peak_amps);
peaks(I(J)) = 0;

peak_indexes = find(peaks);
