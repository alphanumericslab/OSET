function [hr, peaks, hr_period_std] = hr_calculation(x, f, fs, varargin)
% hr_calculation - Calculate heart rate (HR) from an ECG signal.
%
%   [hr, peaks, hr_period_std] = hr_calculation(x, f, fs, peak_det_method, averaging_method, params)
% 
%   Calculates the heart rate (HR) from an ECG signal. The R-peak detector
%   is based on a peak search (peak_det_local_search) on the single
%   or multi-channel signal or its envelope (see modes of operation), and
%   the HR is calculated as the mean, median, or trimmed-mean of the HR
%   over the data window.
%
% Inputs:
%   x: Vector of input ECG data, or matrix in multi-channel mode (channels x samples); see 'peak_det_method'.
%   f: Approximate ECG beat-rate in Hz.
%   fs: Sampling frequency in Hz.
%   peak_det_method: The method used for peak detection.
%       Available options:
%       - 'LOCAL-PEAKS' (default): Detect R-peaks in the original ECG signal.
%       - 'LOCAL-PEAKS-REFINED': Refine R-peak detection using a power envelope approach.
%       - 'POWER-ENVELOPE': Detect R-peaks using the power envelope of the ECG signal.
%       - 'MULTI-LEAD-POWER-ENVELOPE': Detect R-peaks using power envelopes from a multi-lead ECG input.
%   averaging_method: Method to calculate HR ('mean': mean, 'median': median, or 'trmean': trimmed mean).
%   params: A structure containing parameters in different modes of operation.
%       params.wlen: Length of the window used for refined peak detection methods.
%       params.trim: Number of beats to trim from averaging (required only when averaging_method is 'trmean').
%
% Outputs:
%   HR: Calculated heart rate in beats per minute (BPM).
%   peaks: Detected R-peaks.
%   hr_period_std: Standard deviation of the RR-intervals (in samples).
%
% Notes:
% - It is recommended to remove the signal's baseline wander before R-peak detection.
% - The HR is reported in beats per minute (BPM).
%
% Revision History:
%   2009: First release
%   2023: Combined deprecated versions HRCalculation1, HRCalculation2,
%       HRCalculation3, HRCalculation4
%
% Reza Sameni, 2009-2023
% The Open-Source Electrophysiological Toolbox
% https://github.com/alphanumericslab/OSET

% Check if the 'peak_det_method' argument is provided, otherwise use the default 'LOCAL-SEARCH'
if nargin > 3 && ~isempty(varargin{1})
    peak_det_method = varargin{1};
else
    peak_det_method = 'LOCAL-PEAKS';
end

% Check if the 'averaging_method' argument is provided, otherwise use the default 'mean'
if nargin > 4 && ~isempty(varargin{2})
    averaging_method = varargin{2};
else
    averaging_method = 'mean';
end

% Check if the 'params' argument is provided, and extract parameters
if nargin > 5 && ~isempty(varargin{3})
    params = varargin{3};
end

% Detect R-peaks using the specified peak detection method
switch peak_det_method
    case 'LOCAL-PEAKS'
        peaks = peak_det_local_search(x, f/fs);
    case 'LOCAL-PEAKS-REFINED'
        half_wlen = round(params.wlen/2);
        T = length(x);
        x_power_envelope = filtfilt(ones(1, params.wlen), params.wlen, x.^2);
        peaks0 = peak_det_local_search(x_power_envelope, f/fs, 1);
        I0 = find(peaks0);
        I0_len = length(I0);
        peaks = zeros(size(peaks0));
        for i = 1 : I0_len
            index = max(I0(i) - half_wlen, 1) : min(I0(i) + half_wlen, T);
            [~, ii] = max(abs(x(index)));
            peaks(index(1) + ii - 1) = 1;
        end
    case 'POWER-ENVELOPE'
        x_power_envelope = filtfilt(ones(1, params.wlen), params.wlen, x.^2);
        peaks = peak_det_local_search(x_power_envelope, f/fs, 1);
    case 'MULTI-LEAD-POWER-ENVELOPE'
        x_power_envelope = zeros(size(x));
        for k = 1 : size(x_power_envelope, 1)
            x_power_envelope(k, :) = filtfilt(ones(1, params.wlen), params.wlen, x(k, :).^2);
        end
        x_combined_envelopes = sqrt(sum(x_power_envelope, 1));
        peaks = peak_det_local_search(x_combined_envelopes, f/fs, 1);
    otherwise
        error('undefined peak detection method');
end

I = find(peaks);
rr_intervals = diff(I);

% Calculate HR based on the chosen averaging method
switch averaging_method
    case 'mean' % Mean approach
        hr = 60 * fs / mean(rr_intervals);
        hr_period_std = std(rr_intervals);
    case 'median' % Median approach
        hr = 60 * fs / median(rr_intervals);
        hr_period_std = std(rr_intervals);
    case 'trmean' % Trimmed-mean approach
        rr_intervals = sort(rr_intervals);
        hr = 60 * fs / mean(rr_intervals(params.trim + 1:end - params.trim));
        hr_period_std = std(rr_intervals(params.trim + 1:end - params.trim));
    otherwise
        error('undefined averaging method');
end
