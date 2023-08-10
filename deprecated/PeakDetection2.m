function [peaks, peak_indexes] = PeakDetection2(x, fs, varargin)
% PeakDetection has been deprecated. Use
% peak_detection_modified_pan_tompkins instead (with the new ksigma
% parameter), or replace PeakDetection2(data, fs, wlen, fp1, fp2, th, flag), with:
% peak_detection_modified_pan_tompkins(data, fs, wlen, fp1, fp2, th, ksigma, flag)
% Notice the change in the order of ksigma parameter before flag (can be replaced with []).
error('PeakDetection has been deprecated. Use peak_detection_modified_pan_tompkins instead (with the new ksigma parameter), or replace PeakDetection2(data, fs, wlen, fp1, fp2, th, flag), with peak_detection_modified_pan_tompkins(data, fs, wlen, fp1, fp2, th, [], flag). Notice the change in the order of ksigma parameter before flag (can be replaced with [])');
% [peaks, peak_indexes] = peak_detection_modified_pan_tompkins(x, fs, varargin{:});
