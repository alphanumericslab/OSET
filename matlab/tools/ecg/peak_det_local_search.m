function [peaks, peak_indexes] = peak_det_local_search(x, ff, varargin)
%
% peak_det_local_search - R-peak detector based on local max/min search
%
%   [peaks, peak_indexes] = peak_det_local_search(x, f, flag, num_rounds, hr_update_fraction, omit_close_peaks)
%
%   Inputs:
%       x: Vector of input data
%       f: Approximate ECG beat-rate in Hertz, normalized by the sampling frequency
%       flag: Optional. Search for positive (flag=1) or negative (flag=0) peaks.
%           By default, the maximum absolute value of the signal determines the peak sign.
%       num_rounds (optional): Number of iterations to find the R-peaks, up to 3.
%           Default is 1 (run peak detection only once).
%       hr_update_fraction (optional): median HR is multiplied by this
%           fraction in each iteration, if num_rounds > 1. Default
%           is 1.05
%       omit_close_peaks (optional): omit peaks that are too close after main peak 
%           detection (true/1) or not(false/0). Default is false
%
%   Outputs:
%       peaks: Vector of R-peak impulse train
%       peak_indexes: Vector of R-peak indexes
%
%   Notes:
%       - The R-peaks are found from a peak search in windows of length N, where
%         N corresponds to the R-peak period calculated from the given f. R-peaks
%         with periods smaller than N/2 or greater than N are not detected.
%       - It is recommended to remove the signal baseline wander before R-peak detection.
%
%   Revision History:
%       2006: First release
%       2022: Added multi-iteration feature (remains backward compatible)
%       2023: Renamed from deprecated version PeakDetection()
%
%   Reza Sameni, 2006-2023
%   The Open-Source Electrophysiological Toolbox
%   https://github.com/alphanumericslab/OSET

% Set default values for optional arguments
if nargin > 2 && ~isempty(varargin{1})
    flag = varargin{1};
else
    flag = abs(max(x)) > abs(min(x));
end

if nargin > 3 && ~isempty(varargin{2})
    num_rounds = varargin{2};
    if num_rounds > 3
        error('Number of R-peak detections not supported')
    end
else
    num_rounds = 1;
end

if nargin > 4 && ~isempty(varargin{3})
    hr_update_fraction = varargin{3};
else
    hr_update_fraction = 1.05;
end

if nargin > 5 && ~isempty(varargin{4})
    omit_close_peaks = varargin{4};
else
    omit_close_peaks = 0;
end

% Perform peak detection
[peaks, peak_indexes] = peak_det_simple(x, ff, flag, omit_close_peaks);

% Perform additional iterations if specified
if num_rounds > 1
    for k = 1 : num_rounds - 1
        rr_intervals = diff(peak_indexes);
        ff = hr_update_fraction / median(rr_intervals); % Refined heart rate (in Hz) used for R-peak detection
        [peaks, peak_indexes] = peak_det_simple(x, ff, flag, omit_close_peaks);
    end
end