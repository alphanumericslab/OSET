function [HR, peaks, HRPeriodSD] = HRCalculation3(x, f, fs, trim, wlen, averaging_method)
% HRCalculation3 has been deprecated. Use hr_calculation instead. Check the options in hr_calculation
warning('HRCalculation3 has been deprecated. Use hr_calculation instead. Check the options in hr_calculation');
peak_det_method = 'MULTI-LEAD-POWER-ENVELOPE';
params.trim = trim;
params.trim = wlen;
[HR, peaks, HRPeriodSD] = hr_calculation(x, f, fs, peak_det_method, averaging_method, params);
