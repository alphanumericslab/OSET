function [interpeaks, ref_to_out_lag] = IntermediatePeakDetection(x, peaks, wlen, n_peaks, varargin)
% IntermediatePeakDetection has been deprecated. Use intermediate_peak_detection instead.
warning('IntermediatePeakDetection has been deprecated. Use intermediate_peak_detection instead.');
[interpeaks, ref_to_out_lag] = intermediate_peak_detection(x, peaks, wlen, n_peaks, varargin);