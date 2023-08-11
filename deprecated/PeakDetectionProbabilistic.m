function [peaks, peak_indexes, qrs_likelihood] = PeakDetectionProbabilistic(signal, fs, varargin)
% PeakDetectionProbabilistic has been deprecated. Use peak_det_probabilistic instead.
warning('PeakDetectionProbabilistic has been deprecated. Use peak_det_probabilistic instead.');
[peaks, peak_indexes, qrs_likelihood] = peak_det_probabilistic(signal, fs, varargin{:});