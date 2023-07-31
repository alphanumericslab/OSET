clc
clear
close all;

load('testdata.mat');
fs = 1000; % sampling frequency
template_start = 3420; % beginning of reference beat segment
template_stop = 3490; % end of reference beat segment
wlen = round(0.05*fs); % averaging window length
type = 'mean';
PP_diff_th = 50;
PP_diff_wlen = 3;
average_peak_detection_rate = 1.0;
plotflag = 1;

% [PeakToPeak_corrected, PeakToPeak_differences, matched_peaks_indexes] = MFPeakDetector(testdata, template_start, template_stop, type, wlen, PP_diff_wlen, PP_diff_th, average_peak_detection_rate, fs, plotflag);
[PeakToPeak_corrected, PeakToPeak_differences, matched_peaks_indexes, smoothed_matched_output] = MFPeakDetector2(testdata, testdata(template_start:template_stop), type, wlen, PP_diff_wlen, PP_diff_th, average_peak_detection_rate, fs, plotflag);


t = (0:length(testdata)-1)/fs;

figure
plot(t, testdata);
hold on
plot(t,smoothed_matched_output, 'r');
grid;
xlabel('time(s)');
ylabel('Amplitude(mV)');

figure
plot(matched_peaks_indexes(1:end-1)/fs,60*PeakToPeak_corrected/fs);
grid;
xlabel('time(s)');
ylabel('Heart rate(BPM)');
