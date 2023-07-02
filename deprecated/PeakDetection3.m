function [peaks, mn, r] = PeakDetection3(x, fs, h, th, fmax)
% PeakDetection3 has been deprecated. Use peak_detection_matched_filter instead.
    warning('PeakDetection3 has been deprecated. Use peak_detection_matched_filter instead.');
    [peaks, mn, r] = peak_detection_matched_filter(x, fs, h, th, fmax);
end