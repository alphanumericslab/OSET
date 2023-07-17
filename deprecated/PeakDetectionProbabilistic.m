function [peaks, peak_indexes, qrs_likelihood] = PeakDetectionProbabilistic(signal, fs, varargin)
% PeakDetectionProbabilistic has been deprecated. Use peak_detection_probabilistic instead.
warning('PeakDetectionProbabilistic has been deprecated. Use peak_detection_probabilistic instead.');
[peaks, peak_indexes, qrs_likelihood] = peak_detection_probabilistic(signal, fs, varargin{:});