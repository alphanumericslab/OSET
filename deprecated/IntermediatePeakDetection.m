function [interpeaks, ref_to_out_lag] = IntermediatePeakDetection(x, peaks, wlen, n_peaks, varargin)
% IntermediatePeakDetection has been deprecated. Use intermediate_peak_det instead.
warning('IntermediatePeakDetection has been deprecated. Use intermediate_peak_det instead.');
[interpeaks, ref_to_out_lag] = intermediate_peak_det(x, peaks, wlen, n_peaks, varargin);