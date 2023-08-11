function [HR, peaks, HRPeriodSD] = HRCalculation1(x, f, fs, trim, averaging_method)
% HRCalculation1 has been deprecated. Use hr_calculation instead. Check the options in hr_calculation
warning('HRCalculation1 has been deprecated. Use hr_calculation instead. Check the options in hr_calculation');

peak_det_method = 'LOCAL-PEAKS';
params.trim = trim;
[HR, peaks, HRPeriodSD] = hr_calculation(x, f, fs, peak_det_method, averaging_method, params);
