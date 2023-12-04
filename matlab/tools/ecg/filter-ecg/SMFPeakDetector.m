function [PeakToPeak_corrected, PeakToPeak_differences, matched_peaks_indexes, smoothed_matched_output] = SMFPeakDetector(x, matched_template, type,  wlen, PP_diff_wlen, PP_diff_th, average_peak_det_rate, fs)
% ===========================================================
% Matched Filter R-Peak Detector and Smoothed Heart Rate
%
%  inputs:
% 	 x: input vector.
%    matched_template: the reference template
%    type:
%       'mean': a moving average
%       'median': a moving median
% 	 wlen: moving window length of smoothing the energy envelope of the matched filter output (in samples)
% 	 PP_diff_wlen: number of successive beats used for averaging
% 	 PP_diff_th: difference above these number of samples is considered as erroneous and should be corrected
%    average_peak_det_rate: average peak detection rate
%    fs: sampling rate
%  output:
% 	 PeakToPeak_corrected: smoothed matched filter heart rate
%    PeakToPeak_differences: matched filter heart rate without smoothing
%    matched_peaks_indexes: matched filter peak indexes
%    smoothed_matched_output: matched filter energy envelope
%
% Open Source ECG Toolbox, version 1.0, November 2006
% Released under the GNU General Public License
% Copyright (C) 2006  Reza Sameni
% Sharif University of Technology, Tehran, Iran -- LIS-INPG, Grenoble, France
% reza.sameni@gmail.com

% This program is free software; you can redistribute it and/or modify it
% under the terms of the GNU General Public License as published by the
% Free Software Foundation; either version 2 of the License, or (at your
% option) any later version.
% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
% Public License for more details.


% matched filter template selection
% matched_template = x(1,template_stop : -1 : template_start);

% matched filtering
% matched_out = filtfilt(matched_template, mean(matched_template.^2), x);
matched_out = filter(matched_template, mean(matched_template.^2), x);

% matched filter lag compensation
matched_filter_lag = round(length(matched_template)/2);
matched_out = [matched_out(matched_filter_lag:end) zeros(1, matched_filter_lag-1)];

% smoothing the energy envelope of the matched filter output
% smoothed_matched_output = sqrt(filtfilt(ones(wlen,1), wlen, matched_out.^2));
smoothed_matched_output = sqrt(filter(ones(wlen,1), wlen, matched_out.^2));

% secondary matched filter lag compensation
matched_filter_lag2 = round(wlen/2);
smoothed_matched_output = [smoothed_matched_output(matched_filter_lag2:end) zeros(1, matched_filter_lag2-1)];

% peak detection over energy envelope
matched_peaks = PeakDetection(smoothed_matched_output, average_peak_det_rate/fs, 1);
matched_peaks_indexes = find(matched_peaks);

% postprocess the heart rate
PeakToPeak_differences = diff(matched_peaks_indexes);
if length(PeakToPeak_differences) < PP_diff_wlen
    PP_diff_wlen = length(PeakToPeak_differences);
end
PeakToPeak_differences_smoothed = TrimmedFilter(PeakToPeak_differences, type, PP_diff_wlen);

PeakToPeak_gaps = PeakToPeak_differences - PeakToPeak_differences_smoothed;
erroneous_locs = find(abs(PeakToPeak_gaps) > PP_diff_th);
PeakToPeak_corrected = PeakToPeak_differences;
PeakToPeak_corrected(erroneous_locs) = PeakToPeak_differences_smoothed(erroneous_locs);
