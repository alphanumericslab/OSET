function [peaks, peak_indexes, qrs_likelihood] = PeakDetectionProbabilistic(signal, fs, varargin)
% PeakDetectionProbabilistic has been deprecated. Use peak_det_likelihood instead.
warning('PeakDetectionProbabilistic has been deprecated. Use peak_det_likelihood instead.');
[peaks, peak_indexes, qrs_likelihood] = peak_det_likelihood(signal, fs, varargin{:});