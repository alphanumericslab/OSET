function [peaks, peak_indexes, peak_indexes_consensus, qrs_likelihood] = peak_det_likelihood_long_recs(data, fs, varargin)
% PEAK_DET_LIKELIHOOD_LONG_RECS - Block-wise ECG R-peak detector for long records.
%
% This function is an efficient version of peak_det_likelihood,
% designed for long ECG records. It processes the input data in
% segments (blocks) and combines the results to detect R-peaks across the
% entire record.
%
% Usage:
%   [peaks, peak_indexes, peak_indexes_consensus, qrs_likelihood] = peak_det_likelihood_long_recs(data, fs)
%   [peaks, peak_indexes, peak_indexes_consensus, qrs_likelihood] = peak_det_likelihood_long_recs(data, fs, seg_len_time)
%   [peaks, peak_indexes, peak_indexes_consensus, qrs_likelihood] = peak_det_likelihood_long_recs(data, fs, seg_len_time, pad_len_time)
%   [peaks, peak_indexes, peak_indexes_consensus, qrs_likelihood] = peak_det_likelihood_long_recs(data, fs, seg_len_time, pad_len_time, peak_detector_params)
%
% Inputs:
%   - data: The ECG signal data, with each row representing a lead and each
%     column representing a time sample.
%   - fs: The sampling frequency of the ECG data.
%   - seg_len_time (optional): The length of each processing segment in
%     seconds. Default is 10.0 seconds.
%   - pad_len_time (optional): The padding length at the beginning and end
%     of each segment in seconds, for continuity. Default is 1.0 second.
%   - peak_detector_params (optional): Additional parameters for the peak
%     detector passed to peak_det_likelihood(). See the help for
%     peak_det_likelihood, for details
%
% Outputs:
%   - peaks: The detected R-peaks.
%   - peak_indexes: The indices of the detected R-peaks
%   - peak_indexes_consensus: The corrected R-peak indices.
%   - qrs_likelihood: The likelihood of each sample being a QRS complex
%     (between 0 and 1)
%
% Revision History:
%   2023: First release.
%
% Reza Sameni, 2023
% The Open-Source Electrophysiological Toolbox
% https://github.com/alphanumericslab/OSET

% Check for optional input arguments
if nargin > 2 && ~isempty(varargin{1})
    seg_len_time = varargin{1};
else
    seg_len_time = 10.0; % Default segment length in seconds
    disp(['Default seg_len_time = ', num2str(seg_len_time)])
end

if nargin > 3 && ~isempty(varargin{2})
    pad_len_time = varargin{2};
else
    pad_len_time = 1.0; % Default padding length in seconds
    disp(['Default pad_len_time = ', num2str(pad_len_time)])
end

if nargin > 4 && ~isempty(varargin{3})
    peak_detector_params = varargin{3};
else
    peak_detector_params = [];
end

% Calculate segment and padding lengths in samples
seg_len_samples = round(seg_len_time * fs);
pad_len_samples = round(pad_len_time * fs);

% Ensure the segment length does not exceed the data length
if seg_len_samples > size(data, 2)
    seg_len_samples = size(data, 2);
    disp(['Short signal; seg_len_samples truncated to signal length = ', num2str(seg_len_samples)])
end

% Calculate the number of segments
num_seg = ceil(size(data, 2) / seg_len_samples);

% Initialize variables to store results
peak_indexes = [];
peak_indexes_consensus = [];
qrs_likelihood = [];
peaks = [];

% Process each segment
for segment = 1 : num_seg
    first_sample = (segment - 1) * seg_len_samples + 1;
    last_sample = min(size(data, 2), segment * seg_len_samples);

    % Extract the current segment of data
    data_segment = data(:, first_sample : last_sample);

    % Prepare padding data for the peak detector making sure not to exceed
    % the signal end points
    peak_detector_params.left_pad = data(:, max(1, first_sample - pad_len_samples): first_sample - 1);
    if last_sample < size(data, 2)
        peak_detector_params.right_pad = data(:, last_sample + 1: min(last_sample + pad_len_samples, size(data, 2)));
    else
        peak_detector_params.right_pad = [];
    end

    % Call peak_det_likelihood for the current segment
    [peaks_segment, peak_indexes_segment, peak_indexes_consensus_segment, qrs_likelihood_segment] = peak_det_likelihood(data_segment, fs, peak_detector_params);

    % Concatenate results for the current segment with previous segments
    if ~isempty(peak_indexes_segment)
        peak_indexes = cat(2, peak_indexes, peak_indexes_segment + length(peaks));
    end
    if ~isempty(peak_indexes_consensus_segment)
        peak_indexes_consensus = cat(2, peak_indexes_consensus, peak_indexes_consensus_segment + length(peaks));
    end
    qrs_likelihood = cat(2, qrs_likelihood, qrs_likelihood_segment);
    peaks = cat(2, peaks, peaks_segment);
end
