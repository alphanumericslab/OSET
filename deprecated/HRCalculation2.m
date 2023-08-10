function [HR, peaks, HRPeriodSD] = HRCalculation2(x, f, fs, trim, wlen, averaging_method)
% HRCalculation2 has been deprecated. Use hr_calculation instead. Check the options in hr_calculation
warning('HRCalculation2 has been deprecated. Use hr_calculation instead. Check the options in hr_calculation');
peak_detection_method = 'POWER-ENVELOPE';
params.trim = trim;
params.trim = wlen;
[HR, peaks, HRPeriodSD] = hr_calculation(x, f, fs, peak_detection_method, averaging_method, params);
