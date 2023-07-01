% Sample code for two stage matched filter peak detector and rule based
% correction
%
% Fahimeh Jamshidian Tehrani
% Email: jamshidian.t@cse.shirazu.ac.ir, fahimeh.jt@gmaill.com
%
% Created March 2018
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


clc
clear all
close all

load SampleECGData

selectedCH = 2; % selected channel
flgDbg = 1;
% Matched filter peak detector parameters
wide_wlen = round(fs*0.03); % wide moving window length (in samples)
narrow_wlen = round(fs*0.01); % narrow moving window length (in samples)
average_peak_detection_rate = 2.2; %average FECG peak detection rate
PP_diff_wlen = 5; % number of successive beats used for averaging
PP_diff_th = 100; % difference above these number of samples is considered as erroneous and should be corrected
type = 'median';

wide_template_start = 23538;
wide_template_stop = 23599;
narrow_template_start = 23554;
narrow_template_stop = 23578;

if flgDbg
    figure, plot(wide_template_start-fs:wide_template_stop+fs, ECGdata(selectedCH,wide_template_start-fs:wide_template_stop+fs));
    hold on, plot(wide_template_start:wide_template_stop, ECGdata(selectedCH,wide_template_start:wide_template_stop),'r');
    hold on, plot(narrow_template_start:narrow_template_stop, ECGdata(selectedCH,narrow_template_start:narrow_template_stop),'g');
    legend('ECG data', 'wide matched template', 'narrow matched template')
end

% Wide window template
wide_matched_template = ECGdata(selectedCH, wide_template_stop : -1 : wide_template_start);

% Narrow window template
narrow_matched_template = ECGdata(selectedCH, narrow_template_stop : -1 : narrow_template_start);

[multi_matched_peaks_indexes, PeakToPeak_corrected, peak_location_cdi, PeakToPeak_corrected_Rule] = TwoStageSMFPeakDetectorRuleBasedCorrecction(ECGdata(selectedCH,:), wide_matched_template, narrow_matched_template, type,  wide_wlen, narrow_wlen, PP_diff_wlen, PP_diff_th, average_peak_detection_rate, fs, flgDbg);