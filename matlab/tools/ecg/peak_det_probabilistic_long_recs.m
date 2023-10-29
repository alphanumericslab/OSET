function [peaks, peak_indexes, peak_indexes_consensus, qrs_likelihood] = peak_det_probabilistic_long_recs(data, fs, varargin)

if nargin > 2 && ~isempty(varargin{1})
    seg_len_time = varargin{1};
else
    seg_len_time = 10.0;
    disp(['default seg_len_time = ', num2str(seg_len_time)])
end

if nargin > 3 && ~isempty(varargin{2})
    pad_len_time = varargin{2};
else
    pad_len_time = 1.0;
    disp(['default pad_len_time = ', num2str(pad_len_time)])
end

if nargin > 4 && ~isempty(varargin{3})
    peak_detector_params = varargin{3};
else
    peak_detector_params = [];
end

seg_len_samples = round(seg_len_time * fs);
pad_len_samples = round(pad_len_time * fs);
if seg_len_samples > size(data, 2)
    seg_len_samples = size(data, 2);
    disp(['short signal; seg_len_samples truncated to signal length = ', num2str(seg_len_samples)])
end

num_seg = floor(size(data, 2) / seg_len_samples);
peak_indexes = [];
peak_indexes_consensus = [];
qrs_likelihood = [];
peaks = [];
for segment = 1 : num_seg
    first_sample = (segment - 1) * seg_len_samples + 1;
    last_sample = segment * seg_len_samples;

    data_segment = data(:, first_sample : last_sample);


    peak_detector_params.left_pad = data(:, max(1, first_sample - pad_len_samples): first_sample - 1);
    if last_sample < size(data, 2)
        peak_detector_params.right_pad = data(:, last_sample + 1: min(last_sample + pad_len_samples, size(data, 2)));
    else
        peak_detector_params.right_pad = [];
    end
    
    [peaks_segment, peak_indexes_segment, peak_indexes_consensus_segment, qrs_likelihood_segment] = peak_det_probabilistic(data_segment, fs, peak_detector_params);
    peak_indexes = cat(2, peak_indexes, peak_indexes_segment + length(peaks));
    peak_indexes_consensus = cat(2, peak_indexes_consensus, peak_indexes_consensus_segment + length(peaks));
    qrs_likelihood = cat(2, qrs_likelihood, qrs_likelihood_segment);
    peaks = cat(2, peaks, peaks_segment);
end