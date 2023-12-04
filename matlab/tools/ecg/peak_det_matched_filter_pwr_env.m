function [rr_interval_corrected, rr_interval_differences, matched_peaks_indexes, smoothed_matched_output] = peak_det_matched_filter_pwr_env(x, matched_template, type, wlen, PP_diff_wlen, PP_diff_th, average_peak_det_rate, fs, varargin)
% peak_det_matched_filter_pwr_env - matched filter R-Peak detector
% and smoothed heart rate of the matched filter output's power envelope
%
%   [rr_interval_corrected, rr_interval_differences, matched_peaks_indexes, smoothed_matched_output] = peak_det_matched_filter_pwr_env(x, matched_template, type, wlen, PP_diff_wlen, PP_diff_th, average_peak_det_rate, fs, mode)
%
%   This function implements a Matched Filter-based R-Peak detector and
%   smoothed heart rate calculation. It detects R-peaks by convolving the
%   input signal with a matched filter template, then smoothing the energy
%   envelope of the matched filter output, and finally post-processing the
%   detected peaks to correct for errors.
%
% Inputs:
%   x: Input vector (ECG data).
%   matched_template: The reference template for the matched filter.
%   type: Smoothing type ('mean' or 'median').
%   wlen: Moving window length for smoothing the energy envelope of the matched filter output (in samples).
%   PP_diff_wlen: Number of successive beats used for averaging peak-to-peak differences.
%   PP_diff_th: Threshold for detecting erroneous peak-to-peak differences.
%   average_peak_det_rate: Average peak detection rate (peaks per second).
%   fs: Sampling rate in Hz.
%   mode: Optional. Matching mode: 'CAUSAL' (default) or 'NON-CAUSAL'.
%
% Outputs:
%   rr_interval_corrected: Smoothed matched filter heart rate.
%   rr_interval_differences: Matched filter heart rate without smoothing.
%   matched_peaks_indexes: Matched filter peak indexes.
%   smoothed_matched_output: Matched filter energy envelope.
%
% Notes:
% - This function requires the 'trimmed_filter' MEX function from the C++ codes of OSET.
% - The 'mode' parameter determines the causal ('CAUSAL') or non-causal ('NON-CAUSAL') behavior of the matched filter.
%
% Revision History:
%   2006: First release
%   2023: Replaced deprecated version SMFPeakDetector
%
% Reza Sameni, 2006-2023
% The Open-Source Electrophysiological Toolbox
% https://github.com/alphanumericslab/OSET

% Check if the 'mode' argument is provided, otherwise use the default 'CAUSAL'
if nargin > 8 && ~isempty(varargin{1})
    mode = varargin{1};
else
    mode = 'CAUSAL';
end

% Matched filter template selection
% matched_template = x(1, template_stop : -1 : template_start);

% Matched filtering
if strcmp(mode, 'CAUSAL')
    matched_out = filter(matched_template, mean(matched_template.^2), x);
    % Matched filter lag compensation
    matched_filter_lag = round(length(matched_template) / 2);
    matched_out = [matched_out(matched_filter_lag:end) zeros(1, matched_filter_lag - 1)];
elseif strcmp(mode, 'NON-CAUSAL')
    matched_out = filtfilt(matched_template, mean(matched_template.^2), x);
end

% Smoothing the energy envelope of the matched filter output
if strcmp(mode, 'CAUSAL')
    smoothed_matched_output = sqrt(filter(ones(wlen, 1), wlen, matched_out.^2));
    % Secondary matched filter lag compensation
    matched_filter_lag2 = round(wlen / 2);
    smoothed_matched_output = [smoothed_matched_output(matched_filter_lag2:end) zeros(1, matched_filter_lag2 - 1)];
elseif strcmp(mode, 'NON-CAUSAL')
    smoothed_matched_output = sqrt(filtfilt(ones(wlen, 1), wlen, matched_out.^2));
end

% Peak detection over energy envelope
matched_peaks = peak_det_local_search(smoothed_matched_output, average_peak_det_rate / fs, 1);
matched_peaks_indexes = find(matched_peaks);

% Post-process the heart rate
rr_interval_differences = diff(matched_peaks_indexes);
if length(rr_interval_differences) < PP_diff_wlen
    PP_diff_wlen = length(rr_interval_differences);
end
rr_interval_differences_smoothed = trimmed_filter(rr_interval_differences, type, PP_diff_wlen);

rr_interval_gaps = rr_interval_differences - rr_interval_differences_smoothed;
erroneous_locs = find(abs(rr_interval_gaps) > PP_diff_th);
rr_interval_corrected = rr_interval_differences;
rr_interval_corrected(erroneous_locs) = rr_interval_differences_smoothed(erroneous_locs);
