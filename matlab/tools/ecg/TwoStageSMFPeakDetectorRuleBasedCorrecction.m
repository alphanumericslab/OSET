function [multi_matched_peaks_indexes, PeakToPeak_corrected, peak_location_cdi, PeakToPeak_corrected_Rule]=TwoStageSMFPeakDetectorRuleBasedCorrecction(ECG_data, wide_matched_template, narrow_matched_template, type,  wide_wlen, narrow_wlen, PP_diff_wlen, PP_diff_th, average_peak_det_rate, fs, flgDbg)
% ===========================================================
% Two satege (narrow&wide) matched filter peak detector and Rule-Based Correction
%
% Fahimeh Jamshidian Tehrani
% Email: jamshidian.t@cse.shirazu.ac.ir, fahimeh.jt@gmaill.com
%
% Created March 2018
%
% inputs:
%   ECG_data: input ECG record
%   narrow_matched_template: the reference narrow window template
%   wide_matched_template: the reference wide window template
%   type:
%       'mean': a moving average
%       'median': a moving median
%   wide_wlen: wide moving window length (in samples)
%   narrow_wlen: narrow moving window length (in samples)
%   PP_diff_wlen: number of successive beats used for averaging
%   PP_diff_th: difference above these number of samples is considered as erroneous and should be corrected
%   average_peak_det_rate: average FECG peak detection rate
%   fs: sampling rate
%
% output:
%   multi_matched_peaks_indexes: matched filter peak indexes
%   PeakToPeak_corrected: smoothed matched filter heart rate
%   peak_location_cdi: rule based corrected matched filter peak indexes
%   PeakToPeak_corrected_Rule: rule based corrected matched filter heart rate
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


% Call the matched filter R-peak detector and smoothed heart rate function
% for wide and narrow window templates
[PeakToPeak_corrected_HR_all_wide, PeakToPeak_differences_wide, matched_peaks_indexes_wide, smoothed_matched_output_wide] = SMFPeakDetector(ECG_data, wide_matched_template, type,  wide_wlen, PP_diff_wlen, PP_diff_th, average_peak_det_rate, fs);
[PeakToPeak_corrected_HR_all_narrow, PeakToPeak_differences_narrow, matched_peaks_indexes_narrow, smoothed_matched_output_narrow] = SMFPeakDetector(ECG_data, narrow_matched_template, type,  narrow_wlen, PP_diff_wlen, PP_diff_th, average_peak_det_rate, fs);

% peak detection over energy envelope
multi_matched_peaks = PeakDetection(smoothed_matched_output_wide.*smoothed_matched_output_narrow, average_peak_det_rate/fs, 1);
multi_matched_peaks_indexes = find(multi_matched_peaks);

% postprocess the heart rate 
multi_PeakToPeak_differences = diff(multi_matched_peaks_indexes);
PeakToPeak_differences_smoothed = TrimmedFilter(multi_PeakToPeak_differences, type, PP_diff_wlen);
PeakToPeak_gaps = multi_PeakToPeak_differences - PeakToPeak_differences_smoothed;
erroneous_locs = find(abs(PeakToPeak_gaps) > PP_diff_th);
PeakToPeak_corrected = multi_PeakToPeak_differences;
PeakToPeak_corrected(erroneous_locs) = PeakToPeak_differences_smoothed(erroneous_locs);

% Rule based R-peak correction
[peak_location_c, peak_location_cd, peak_location_cdi] = correctPeakLocation(multi_matched_peaks_indexes, ECG_data, flgDbg);
PeakToPeak_corrected_Rule = diff( peak_location_cdi);
end

function  [peak_location_c, peak_location_cd, peak_location_cdi] = correctPeakLocation(peak_location, ECG_data, flgDbg)
% ===========================================================
% The rule based peak location correction 
% The correction is based on heart rate variation
%
% input:
% peak_location: detected peak location 
% ECG_data: data for plot in debug mode
% output:
% peak_location_c: corrected peak location in which 
%                  the wrong located R-peaks are corrected
% peak_location_cd: corrected peak location in which 
%                   the wrong located R-peaks are corrected and
%                   the extra detected R-peaks are deleted
% peak_location_cdi: corrected peak location in which 
%                   the wrong located R-peaks are corrected and
%                   the extra detected R-peaks are deleted and
%                   the missing R-peaks are inserted 
%  


CorrectedPeakLocation = peak_location;
T = mean(diff(peak_location));
% Correcting the wrong location R-peak(s)
PeakToPeak_2differences = diff(CorrectedPeakLocation, 2);
pow_PeakToPeak_2differences = PeakToPeak_2differences(1:end-1) .* abs(PeakToPeak_2differences(2:end));
pow_pow = pow_PeakToPeak_2differences(1:end-1) .* pow_PeakToPeak_2differences(2:end);
normal_pow_pow = (pow_pow).*std(PeakToPeak_2differences)./std(pow_pow);
errorInd = find(normal_pow_pow < -T/8);
CorrectedPeakLocation(errorInd+2) = round(mean([CorrectedPeakLocation(errorInd+1)', CorrectedPeakLocation(errorInd+3)'],2)');
peak_location_c = CorrectedPeakLocation;

% Deleting extra detected R-peak(s)
PeakToPeak_2differences = diff(CorrectedPeakLocation,2);
ind_2differences = find(PeakToPeak_2differences<-T/5)+2;
PeakToPeak_differences = (diff(CorrectedPeakLocation));
T = mean(PeakToPeak_differences);
th = T-.2*T;
ind_differences = find(PeakToPeak_differences<th);
dell_ind = ind_2differences(ismember(ind_2differences,ind_differences));
CorrectedPeakLocation(dell_ind) = [];
peak_location_cd = CorrectedPeakLocation;

%Insert missing R-peak(s) and delete the wrong detected R-peak(s)
PeakToPeak_2differences = diff(CorrectedPeakLocation,2);
ind_2differences = find(PeakToPeak_2differences>T/2.5)+2;
PeakToPeak_differences = diff(CorrectedPeakLocation);
th = T+.2*T;
ind_differences = find(PeakToPeak_differences>th);
cpltemp = CorrectedPeakLocation;
insert_ind = ind_2differences(ismember(ind_2differences,ind_differences));
dif = cpltemp(insert_ind+1)-cpltemp(insert_ind-1);
firstPeak =  round(cpltemp(insert_ind-1) + dif*1/3);
secondPeak =  round(cpltemp(insert_ind-1) + dif*2/3);
CorrectedPeakLocation(insert_ind) = [];
CorrectedPeakLocation = sort([CorrectedPeakLocation firstPeak secondPeak]);
peak_location_cdi = CorrectedPeakLocation;

if flgDbg
    figure,
    plot(diff(peak_location))
    hold on
    plot(diff(peak_location_c))
    plot(diff(peak_location_cd))
    plot(diff(peak_location_cdi))
    legend('heart rate(HR)','HR after correcting wrong location','HR after deleting extra R-peaks','HR after inserting missing R-peaks')
    
    figure,
    plot(ECG_data)
    hold on
    plot(peak_location, ECG_data(peak_location),'ro')
    plot(peak_location_c, ECG_data(peak_location_c),'gs')
    plot(peak_location_cd, ECG_data(peak_location_cd),'m*')
    plot(peak_location_cdi, ECG_data(peak_location_cdi),'kd')
    legend('ECG data','peak location(PL)','PL after correcting wrong location','PL after deleting extra R-peaks','PL after inserting missing R-peaks')
end
end


