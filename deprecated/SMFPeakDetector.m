function [PeakToPeak_corrected, PeakToPeak_differences, matched_peaks_indexes, smoothed_matched_output] = SMFPeakDetector(x, matched_template, type,  wlen, PP_diff_wlen, PP_diff_th, average_peak_det_rate, fs)
% SMFPeakDetector has been deprecated. Use peak_det_matched_filter_pwr_env instead.
warning('SMFPeakDetector has been deprecated. Use peak_det_matched_filter_pwr_env instead.');
[PeakToPeak_corrected, PeakToPeak_differences, matched_peaks_indexes, smoothed_matched_output] = peak_det_matched_filter_pwr_env(x, matched_template, type,  wlen, PP_diff_wlen, PP_diff_th, average_peak_det_rate, fs);
